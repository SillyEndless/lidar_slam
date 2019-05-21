#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <iostream>


//位姿-->转换矩阵
Eigen::Matrix3d PoseToTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d trans;
    trans << cos(x(2)),-sin(x(2)),x(0),
             sin(x(2)), cos(x(2)),x(1),
                     0,         0,    1;

    return trans;
}


//转换矩阵－－＞位姿
Eigen::Vector3d TransToPose(Eigen::Matrix3d trans)
{
    Eigen::Vector3d pose;
    pose(0) = trans(0,2);
    pose(1) = trans(1,2);
    pose(2) = atan2(trans(1,0),trans(0,0));

    return pose;
}

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<Edge>& Edges)
{
    double sumError = 0;
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z  = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() *  Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);


        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}


/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d xi,Eigen::Vector3d xj,Eigen::Vector3d z,
                          Eigen::Vector3d& ei,Eigen::Matrix3d& Ai,Eigen::Matrix3d& Bi)
{

    Eigen::Matrix3d Xi = PoseToTrans(xi);
    Eigen::Matrix3d Xj = PoseToTrans(xj);
    Eigen::Matrix3d Zij  = PoseToTrans(z);

    Eigen::Matrix3d Eij = Zij.inverse() *  Xi.inverse() * Xj;

    ei = TransToPose(Eij);//观测值和预测值之间的误差
    //以下为求解Ai,Bi所需参数。
    Eigen::Matrix2d Rij = Zij.block(0,0,2,2);
    Eigen::Matrix2d Ri = Xi.block(0,0,2,2);

    Eigen::Vector2d ti,tj;
    ti << xi(0),xi(1);
    tj << xj(0),xj(1);

    double theta = xi(2);

    Eigen::Matrix2d dRiT;//Ri对theta求导后的转置
    dRiT << -sin(theta),cos(theta),
            -cos(theta),-sin(theta);

    //以下为使用公式
    Ai.setZero();
    Ai.block(0,0,2,2) = - Rij.transpose() * Ri.transpose();
    Ai.block(0,2,2,1) = Rij.transpose() * dRiT * (tj - ti);
    Ai(2,2) = -1;
    Bi.setZero();
    Bi.block(0,0,2,2) = Rij.transpose() * Ri.transpose();
    Bi(2,2) = 1;

}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd  LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<Edge>& Edges)
{
    //申请内存
    Eigen::MatrixXd H(Vertexs.size() * 3,Vertexs.size() * 3);
    Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0,0,3,3) += I;/*即H11为 1,0,0
                                   0,1,0,
                                   0,0,1
                                  确保dx1为0,0,0 */
    //构造H矩阵　＆ b向量
    for(int k = 0; k < Edges.size();k++)
    {


        //提取信息
        Edge tmpEdge = Edges[k];//取第k条边
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];//得到连接这条边的两个节点（的位姿）
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d zij = tmpEdge.measurement;//得到这两个节点之间的观测值。节点j相对于节点i的位姿。帧间匹配得到。还有回环。
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;//匹配时候得到，协方差矩阵的逆。

        //计算误差和对应的Jacobian
        Eigen::Vector3d eij;
        Eigen::Matrix3d Aij;
        Eigen::Matrix3d Bij;
        CalcJacobianAndError(xi,xj,zij,eij,Aij,Bij);

        Eigen::Matrix3d Hii,Hij,Hji,Hjj;
        Hii = Aij.transpose() * infoMatrix * Aij;
        Hij = Aij.transpose() * infoMatrix * Bij;
        Hji = Bij.transpose() * infoMatrix * Aij;
        Hjj = Bij.transpose() * infoMatrix * Bij;

        Eigen::Vector3d bi,bj;
        bi = (eij.transpose() * infoMatrix * Aij).transpose();
        bj = (eij.transpose() * infoMatrix * Bij).transpose();

        int ix = tmpEdge.xi;
        int jx = tmpEdge.xj;
        //H11为单位矩阵。给的数据有第0个节点，这里表示为H00为单位阵更好理解。
        H.block(3*ix,3*ix,3,3) += Hii;
        H.block(3*ix,3*jx,3,3) += Hij;
        H.block(3*jx,3*ix,3,3) += Hji;
        H.block(3*jx,3*jx,3,3) += Hjj;
        //b的前3维为0,0,0,
        b(3*ix) += bi(0);
        b(3*ix+1) += bi(1);
        b(3*ix+2) += bi(2);
        b(3*jx) += bj(0);
        b(3*jx+1) += bj(1);
        b(3*jx+2) += bj(2);

    }

    //求解
    Eigen::VectorXd dx;


    //dx是3×n维向量，n是图优化中节点的个数。
    //每三维是一个相对位姿优化。
    //dx1肯定为0,0,0.(由H11为单位矩阵保证.)
    dx = -H.lu().solve(b);

    return dx;

}











