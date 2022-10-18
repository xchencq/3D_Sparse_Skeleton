#include "polygon_generation/skeleton_finder_3D.h"

int main(int argc, char **argv) {
  ros::init(argc, argv, "skeleton_finder_3D_test");
  ros::NodeHandle n("~");

  SkeletonFinder skeleton_finder_3D;
  skeleton_finder_3D.init(n);

  ros::Duration(1.0).sleep();
  ros::spin();
}