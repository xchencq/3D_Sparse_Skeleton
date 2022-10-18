#include <math.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <sys/time.h>
#include <time.h>
#include <pcl/search/impl/kdtree.hpp>

using namespace std;
string file_name;

int main(int argc, char **argv) {
  ros::init(argc, argv, "map generator");
  ros::NodeHandle nodehandle("~");
  ros::Rate loop_rate(1);

  file_name = argv[1];

  ros::Publisher map_pub = nodehandle.advertise<sensor_msgs::PointCloud2>("map", 1);

  sensor_msgs::PointCloud2 map_pcd;

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);

  // Read map
  pcl::PCDReader reader;
  reader.read(file_name, *cloud);

  ROS_INFO("Size of map: %d", (*cloud).points.size());
  pcl::toROSMsg(*cloud, map_pcd);
  map_pcd.header.frame_id = "map";

  while (ros::ok()) {
    map_pub.publish(map_pcd);
    ros::spinOnce();
    loop_rate.sleep();
  }
}