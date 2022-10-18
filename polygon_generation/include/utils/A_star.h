#include <math.h>
#include <Eigen/Eigen>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <queue>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/search/impl/kdtree.hpp>

#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include "polygon_generation/backward.hpp"
#include "polygon_generation/data_type_3D.h"

namespace a_star {
struct Node;
typedef Node* NodePtr;

struct Node {
  Eigen::Vector3d coord;
  double g_score, f_score;
  bool expanded;
  NodePtr parent;
  vector<NodePtr> connected_node;
  vector<Eigen::Vector3d> connected_node_pos;

  Node(Eigen::Vector3d coord_, vector<Eigen::Vector3d> connected_node_pos_) {
    coord = coord_;
    connected_node_pos = connected_node_pos_;

    // Set default values;
    g_score = -1;
    f_score = -1;
    expanded = false;
    parent = NULL;
  }
};

struct NodeComparator {
public:
  bool operator()(NodePtr node1, NodePtr node2) {
    return node1->f_score > node2->f_score;
  }
};

class AStar {
 public:
  // Parameter
  double lambda_heu_;

  AStar();
  ~AStar();
  enum { REACH_END = 1, NO_PATH = 2 };

  void initHandler(ros::NodeHandle &n);
  void init(vector<NodePtr> nodes);
  static double pathLength(const vector<Eigen::Vector3d>& path);
  int search(const Eigen::Vector3d& start_pt, const Eigen::Vector3d& end_pt);

  std::vector<Eigen::Vector3d> getPath();
  void visNodes();
  void visConnections();

 private:
  // ros related
  ros::Publisher vis_node_pub, vis_connections_pub, vis_parents_pub;
  pcl::PointCloud<pcl::PointXYZ> nodes_pcl;
  sensor_msgs::PointCloud2 nodes_pcl_ros;

  // parameter
  double tie_breaker_;

  // main data structure
  vector<NodePtr> node_pool;
  std::priority_queue<NodePtr, std::vector<NodePtr>, NodeComparator> queue;
  std::vector<Eigen::Vector3d> path_nodes;

  void backtrack(const NodePtr& end_node, const Eigen::Vector3d& end);

  double getEuclHeu(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2);
  double getEuclDis(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2);
  NodePtr getNodeFromPos(Eigen::Vector3d pos);
  bool isSamePos(Eigen::Vector3d p1, Eigen::Vector3d p2);
  bool isEndNode(NodePtr node, Eigen::Vector3d end);
  void visParents();

  
};
}  // namespace a_star