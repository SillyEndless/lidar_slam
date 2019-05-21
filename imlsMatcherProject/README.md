1.首先编译安装 champion_nav_msgs,按照 champion_nav_msgs 的 readme 文件执行即可,注意根据自己的ros版本进行修改。
2.在imlsMatcherProject目录下catkin_make
3.source devel/setup.bash
4.roscore
5.rosrun imlsMatcher imlsMatcher_node 
6.下载数据集.
链接：https://pan.baidu.com/s/1ECdub5m7JHwPLNvGxH3Vhw 
提取码：xczt 
并执行rosbag play --clock imls_icp.bag

