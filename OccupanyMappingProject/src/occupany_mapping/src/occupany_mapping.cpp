#include "occupany_mapping.h"
#include "nav_msgs/GetMap.h"
#include "sensor_msgs/PointCloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "geometry_msgs/Point32.h"


/**
 * Increments all the grid cells from (x0, y0) to (x1, y1);
 * //不包含(x1,y1)
 * 2D画线算法　来进行计算两个点之间的grid cell
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 */
std::vector<GridIndex> TraceLine(int x0, int y0, int x1, int y1)
{
  GridIndex tmpIndex;
  std::vector<GridIndex> gridIndexVector;

  bool steep = abs(y1 - y0) > abs(x1 - x0);
  if (steep)
  {
    std::swap(x0, y0);
    std::swap(x1, y1);
  }
  if (x0 > x1)
  {
    std::swap(x0, x1);
    std::swap(y0, y1);
  }

  int deltaX = x1 - x0;
  int deltaY = abs(y1 - y0);
  int error = 0;
  int ystep;
  int y = y0;

  if (y0 < y1)
  {
    ystep = 1;
  }
  else
  {
    ystep = -1;
  }

  int pointX;
  int pointY;
  for (int x = x0; x <= x1; x++)
  {
    if (steep)
    {
      pointX = y;
      pointY = x;
    }
    else
    {
      pointX = x;
      pointY = y;
    }

    error += deltaY;

    if (2 * error >= deltaX)
    {
      y += ystep;
      error -= deltaX;
    }

    //不包含最后一个点．
    if(pointX == x1 && pointY == y1) continue;

    //保存所有的点
    tmpIndex.SetIndex(pointX,pointY);

    gridIndexVector.push_back(tmpIndex);
  }

  return gridIndexVector;
}

void SetMapParams(void )
{
   mapParams.width = 900;
   mapParams.height = 900;
   mapParams.resolution = 0.04;

   //每次被击中的log变化值
   mapParams.log_free = -1;
   mapParams.log_occ = 2;

   //每个栅格的最大最小值．
   mapParams.log_max = 100.0;
   mapParams.log_min = 0.0;

   mapParams.origin_x = 0.0;
   mapParams.origin_y = 0.0;

   //地图的原点，在地图的正中间
   mapParams.offset_x = 700;
   mapParams.offset_y = 600;

   pMap = new unsigned char[mapParams.width*mapParams.height];

   //初始化为50
   for(int i = 0; i < mapParams.width * mapParams.height;i++)
        pMap[i] = 50;
}


//从世界坐标系转换到栅格坐标系
GridIndex ConvertWorld2GridIndex(double x,double y)
{
    GridIndex index;

    index.x = std::ceil((x - mapParams.origin_x) / mapParams.resolution) + mapParams.offset_x;
    index.y = std::ceil((y - mapParams.origin_y) / mapParams.resolution) + mapParams.offset_y;

    return index;
}

int GridIndexToLinearIndex(GridIndex index)
{
    int linear_index;
    linear_index = index.y + index.x * mapParams.width;
    return linear_index;
}


//判断index是否有效
bool isValidGridIndex(GridIndex index)
{
    if(index.x >= 0 && index.x < mapParams.width && index.y >= 0 && index.y < mapParams.height)
        return true;

    return false;
}

void DestoryMap()
{
    if(pMap != NULL)
        delete pMap;
}

//
void OccupanyMapping(std::vector<GeneralLaserScan>& scans,std::vector<Eigen::Vector3d>& robot_poses)
{
    //枚举所有的激光雷达数据
    for(int i = 0; i < scans.size();i++)
    {
        GeneralLaserScan scan = scans[i];
        Eigen::Vector3d robotPose = robot_poses[i];

        //机器人的地图坐标系
        GridIndex robotIndex = ConvertWorld2GridIndex(robotPose(0),robotPose(1));

        for(int id = 0; id < scan.range_readings.size();id++)
        {
            double dist = scan.range_readings[id];
            double angle = -scan.angle_readings[id]; // 激光雷达逆时针转，角度取反

            if(std::isinf(dist) || std::isnan(dist)) continue;

            //计算得到第i帧激光雷达数据中的第id束激光点的世界坐标系的坐标
            double theta = -robotPose(2); // 激光雷达逆时针转，角度取反
            double laser_x = dist * cos(angle);
            double laser_y = dist * sin(angle);

            double world_x = cos(theta) * laser_x - sin(theta) * laser_y + robotPose(0);
            double world_y = sin(theta) * laser_x + cos(theta) * laser_y + robotPose(1);

            //对对应的map的cell进行更新．
            GridIndex GIndex = ConvertWorld2GridIndex(world_x,world_y);//转换到地图坐标系

            if(isValidGridIndex(GIndex) == false||isValidGridIndex(robotIndex) == false)continue;

            //得到机器人与激光点之间所有的被该束激光通过的栅格的坐标
            std::vector<GridIndex> freeIndexs = TraceLine(robotIndex.x,robotIndex.y,GIndex.x,GIndex.y);
            //更新这些被通过的栅格
            for(int k = 0; k < freeIndexs.size();k++)
            {
                GridIndex tmpIndex = freeIndexs[k];
                int linearIndex = GridIndexToLinearIndex(tmpIndex);

                int tmpData = pMap[linearIndex];
                tmpData += mapParams.log_free;
                if(tmpData < 0)tmpData = 0;
                pMap[linearIndex] = tmpData;
            }
            //更新被击中的点．
            int linearIndex = GridIndexToLinearIndex(GIndex);//激光点的地图坐标。激光点就是激光束击中的障碍物
            int tmpData = pMap[linearIndex];
            tmpData += mapParams.log_occ;
            if(tmpData > 100)tmpData = 100;
            pMap[linearIndex] = tmpData;

        }
    }
}


//发布地图．
void PublishMap(ros::Publisher& map_pub)
{
    nav_msgs::OccupancyGrid rosMap;

    rosMap.info.resolution = mapParams.resolution;
    rosMap.info.origin.position.x = 0.0;
    rosMap.info.origin.position.y = 0.0;
    rosMap.info.origin.position.z = 0.0;
    rosMap.info.origin.orientation.x = 0.0;
    rosMap.info.origin.orientation.y = 0.0;
    rosMap.info.origin.orientation.z = 0.0;
    rosMap.info.origin.orientation.w = 1.0;

    rosMap.info.origin.position.x = mapParams.origin_x;
    rosMap.info.origin.position.y = mapParams.origin_y;
    rosMap.info.width = mapParams.width;
    rosMap.info.height = mapParams.height;
    rosMap.data.resize(rosMap.info.width * rosMap.info.height);

    //0~100
    int cnt0,cnt1,cnt2;
    cnt0 = cnt1 = cnt2 = 100;
    for(int i = 0; i < mapParams.width * mapParams.height;i++)
    {
       if(pMap[i] == 50)
       {
           rosMap.data[i] = -1.0;
       }
       else
       {

           rosMap.data[i] = pMap[i];
       }
    }

    rosMap.header.stamp = ros::Time::now();
    rosMap.header.frame_id = "map";

    map_pub.publish(rosMap);
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "OccupanyMapping");

    ros::NodeHandle nodeHandler;

    ros::Publisher mapPub = nodeHandler.advertise<nav_msgs::OccupancyGrid>("laser_map",1,true);

    std::vector<Eigen::Vector3d> robotPoses;
    std::vector<GeneralLaserScan> generalLaserScans;

    std::string basePath = "/home/elliott/Documents/Lslam/OccupanyMappingProject/src/occupany_mapping/data";

    std::string posePath= basePath + "/pose.txt";
    std::string anglePath = basePath + "/scanAngles.txt";
    std::string scanPath = basePath + "/ranges.txt";

    //读取数据
    ReadPoseInformation(posePath,robotPoses);

    ReadLaserScanInformation(anglePath,
                             scanPath,
                             generalLaserScans);

    //设置地图信息
    SetMapParams();

    OccupanyMapping(generalLaserScans,robotPoses);

    PublishMap(mapPub);

    ros::spin();

    DestoryMap();

    std::cout <<"Release Memory!!"<<std::endl;

}




