#include <map.h>
#include "gaussian_newton_method.h"
#include <math.h>       /* sqrt */
const double GN_PI = 3.1415926;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
          sin(vec(2)), cos(vec(2)),vec(1),
                    0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 10000;
    map->size_y = 10000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);//将激光点在激光参考系下的坐标转换为在里程计下的坐标。

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);//障碍物影响的最大距离为0.5

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示市场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;

    double tmp_x,tmp_y;
    tmp_x=coords[0];
    tmp_y=coords[1];//里程计坐标（世界坐标）
    double x = MAP_GXWX(map,tmp_x);
    double y = MAP_GYWY(map,tmp_y);//栅格坐标

    int x0,x1,y0,y1;
    x0 = x;x1 = x + 1;y0 = y;y1 = y+1;

    x=(tmp_x - map->origin_x) / map->resolution + 0.5 + map->size_x / 2;
    y=(tmp_y - map->origin_y) / map->resolution + 0.5 + map->size_y / 2;//地图坐标
    //检测栅格坐标是否出界
    if (!MAP_VALID(map, x0, y0)||!MAP_VALID(map, x1, y0)
        ||!MAP_VALID(map, x0, y1)||!MAP_VALID(map, x1, y1))
        return Eigen::Vector3d(0.0d,0.0d,0.0d);
    double MP[2][2];
    MP[0][0] = map->cells[MAP_INDEX(map,x0,y0)].score;
    MP[1][0] = map->cells[MAP_INDEX(map,x1,y0)].score;
    MP[0][1] = map->cells[MAP_INDEX(map,x0,y1)].score;
    MP[1][1] = map->cells[MAP_INDEX(map,x1,y1)].score;

    ans[0]=(y-y0)*((x-x0)*MP[1][1]+(x1-x)*MP[0][1])+
           (y1-y)*((x-x0)*MP[1][0]+(x1-x)*MP[0][0]);
    ans[1]=(y-y0)*(MP[1][1]-MP[0][1])+(y1-y)*(MP[1][0]-MP[0][0]);
    ans[2]=(x-y0)*(MP[1][1]-MP[1][0])+(x1-x)*(MP[0][1]-MP[0][0]);
    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
void ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                                  std::vector<Eigen::Vector2d>& laser_pts,
                                  Eigen::Matrix3d& H, Eigen::Vector3d& dTr)
{
    H = Eigen::Matrix3d::Zero();
    dTr = Eigen::Vector3d::Zero();

    int size = laser_pts.size();

    now_pose(2) = GN_NormalizationAngle(now_pose(2));
    double sinRot = sin(now_pose[2]);//now_pose[2]是Ttheta
    double cosRot = cos(now_pose[2]);

    Eigen::Matrix3d Trans = GN_V2T(now_pose);//now_pose是激光发射器在里程计坐标系下的位姿。

    for (int i = 0; i < size; ++i) {

        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);//将激光点在当前激光参考系下的坐标转换为在里程计下的坐标。
        Eigen::Vector3d transformedPointData(InterpMapValueWithDerivatives(map,tmp_pt));//map是上一帧激光构建的栅格地图。

        double funVal = 1.0d - transformedPointData[0];
/* dTr 这个向量计算的是公式12 的求和符号后面的部分，没有乘H-1*/
        dTr[0] += transformedPointData[1] * funVal;
        dTr[1] += transformedPointData[2] * funVal;
/*
这个对应论文中公式 13，14 ,计算 M的梯度叉乘S对旋转角度的偏导数
*/
//        double rotDeriv = ((-sinRot * tmp_pt.x() - cosRot * tmp_pt.y())
//                          * transformedPointData[1] +
//                          (cosRot * tmp_pt.x() - sinRot * tmp_pt.y())
//                          * transformedPointData[2]);
        double rotDeriv = ((-sinRot * laser_pts[i].x() - cosRot * laser_pts[i].y())
                           * transformedPointData[1] +
                           (cosRot * laser_pts[i].x() - sinRot * laser_pts[i].y())
                           * transformedPointData[2]);//这里使用laser_pts[i].x()，激光点相对于激光发射器的坐标。

        dTr[2] += rotDeriv * funVal;

        H(0, 0) += pow(transformedPointData[1],2);
        H(1, 1) += pow(transformedPointData[2],2);
        H(2, 2) += pow(rotDeriv,2);
        H(0, 1) += transformedPointData[1] * transformedPointData[2];
        H(0, 2) += transformedPointData[1] * rotDeriv;
        H(1, 2) += transformedPointData[2] * rotDeriv;
    }

    H(1, 0) = H(0, 1);
    H(2, 0) = H(0, 2);
    H(2, 1) = H(1, 2);
}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map map是由上一帧激光数据构建的地图
 * @param init_pose 初始值。
 * @param laser_pts
 *
 * 利用初始值，不断进行迭代优化，20次后。
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    int maxIteration = 20;
    Eigen::Vector3d now_pose = init_pose;//
    Eigen::Matrix3d H;
    Eigen::Vector3d b;
    Eigen::Vector3d dT;
    for(int i = 0; i < maxIteration;i++)
    {
        now_pose(2) = GN_NormalizationAngle(now_pose(2));
        ComputeHessianAndb(map,now_pose,laser_pts,H,b);
        if ((H(0, 0) == 0.0f) || (H(1, 1) == 0.0f)) continue;
        dT=H.inverse()*b;
        dT(2) = GN_NormalizationAngle(dT(2));
        now_pose+=dT*0.1;

    }
    init_pose = now_pose;

    std::cout <<"GaussianNewton--end_pose:"<<init_pose(0)<<","<<init_pose(1)<<","<<init_pose(2)*57.295<<std::endl;
}

