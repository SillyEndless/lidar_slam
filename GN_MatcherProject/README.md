1.在GN_MatcherProject目录下catkin_make
2.source devel/setup.bash
3.roscore
4.rosrun gaussian_newton_scanmatcher gaussian_newton_node 

5.下载数据集.
链接：https://pan.baidu.com/s/1djUoPS4HuBVsHcNaB_xK0Q 
提取码：68ff 
并执行rosbag play --clock odom.bag

6.启动rviz。观察矫正前后路径。(根据代码中topic配置rviz)

