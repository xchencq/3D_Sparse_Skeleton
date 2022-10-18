#include "utils/A_star.h"
#include <algorithm>

// using namespace Eigen;
// using namespace std;

// namespace backward {
// backward::SignalHandling sh;
// }

namespace a_star {
AStar::AStar() {}

AStar::~AStar() {}

double AStar::pathLength(const vector<Eigen::Vector3d>& path) {
  double length = 0.0;
  if (path.size() < 2) return length;
  for (int i = 0; i < path.size() - 1; ++i) length += (path[i + 1] - path[i]).norm();
  return length;
}

void AStar::initHandler(ros::NodeHandle& n) {
  vis_node_pub = n.advertise<sensor_msgs::PointCloud2>("/A_star/nodes", 1);
  vis_connections_pub = n.advertise<visualization_msgs::Marker>("/A_star/connections", 1);
  vis_parents_pub = n.advertise<visualization_msgs::Marker>("/A_star/parents", 1);
}

void AStar::init(vector<NodePtr> nodes) {
  tie_breaker_ = 1.0 + 1.0 / 1000;
  for (NodePtr node : nodes) {
    node_pool.push_back(node);
  }

  for (NodePtr node : nodes) {
    for (Eigen::Vector3d con_node_pos : node->connected_node_pos) {
      NodePtr con_node = getNodeFromPos(con_node_pos);
      if (con_node == NULL) {
        ROS_ERROR("Can't find Node at given position(%f, %f, %f)!", con_node_pos(0),
                  con_node_pos(1), con_node_pos(2));
        continue;
      }
      node->connected_node.push_back(con_node);
    }
  }
}

int AStar::search(const Eigen::Vector3d& start_pt, const Eigen::Vector3d& end_pt) {
  NodePtr start = getNodeFromPos(start_pt);
  start->g_score = 0;
  start->f_score = getEuclHeu(start_pt, end_pt);
  queue.push(start);

  NodePtr cur_node;
  while (!queue.empty()) {
    cur_node = queue.top();
    // ROS_WARN("Expanding node (%f, %f, %f), g = %f, f = %f", cur_node->coord(0), cur_node->coord(1),
    //          cur_node->coord(2), cur_node->g_score, cur_node->f_score);
    if (isEndNode(cur_node, end_pt)) {
      backtrack(cur_node, end_pt);
      return REACH_END;
    }

    queue.pop();
    cur_node->expanded = true;

    // ROS_INFO("Number of connected nodes: %d", cur_node->connected_node.size());
    for (NodePtr nbhd_node : cur_node->connected_node) {
      if (nbhd_node->expanded) continue;

      double nbhd_g = cur_node->g_score + getEuclDis(cur_node->coord, nbhd_node->coord);
      if (nbhd_node->g_score < 0) {
        nbhd_node->parent = cur_node;
        // ROS_INFO("Change parent of node (%f, %f, %f) to current node.", nbhd_node->coord(0),
        //          nbhd_node->coord(1), nbhd_node->coord(2));
        nbhd_node->g_score = nbhd_g;
        nbhd_node->f_score = nbhd_node->g_score + getEuclHeu(nbhd_node->coord, end_pt);
        queue.push(nbhd_node);
      } else if (nbhd_node->g_score > nbhd_g) {
        nbhd_node->g_score = nbhd_g;
        nbhd_node->parent = cur_node;
        // ROS_INFO("Change parent of node (%f, %f, %f) to current node.", nbhd_node->coord(0),
        //          nbhd_node->coord(1), nbhd_node->coord(2));
      }
    }

    // ROS_WARN("-------------------------------------");
    // visParents();
    // getchar();
  }

  return NO_PATH;
}

std::vector<Eigen::Vector3d> AStar::getPath() { return path_nodes; }

void AStar::backtrack(const NodePtr& end_node, const Eigen::Vector3d& end) {
  path_nodes.push_back(end);
  path_nodes.push_back(end_node->coord);
  NodePtr cur_node = end_node;
  while (cur_node->parent != NULL) {
    cur_node = cur_node->parent;
    path_nodes.push_back(cur_node->coord);
  }
  reverse(path_nodes.begin(), path_nodes.end());
}

double AStar::getEuclHeu(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) {
  return tie_breaker_ * (x2 - x1).norm();
}

double AStar::getEuclDis(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) {
  return (x2 - x1).norm();
}

NodePtr AStar::getNodeFromPos(Eigen::Vector3d pos) {
  for (NodePtr node : node_pool) {
    if (isSamePos(node->coord, pos)) return node;
  }
  return NULL;
}

bool AStar::isSamePos(Eigen::Vector3d p1, Eigen::Vector3d p2) {
  return fabs(p1(0) - p2(0)) < 1e-3 && fabs(p1(1) - p2(1)) < 1e-3 && fabs(p1(2) - p2(2)) < 1e-3;
}

bool AStar::isEndNode(NodePtr node, Eigen::Vector3d end) { return isSamePos(node->coord, end); }

void AStar::visNodes() {
  nodes_pcl.clear();
  int num_nodes = node_pool.size();
  for (int i = 0; i < num_nodes; i++) {
    NodePtr node = node_pool.at(i);
    nodes_pcl.points.push_back(pcl::PointXYZ(node->coord(0), node->coord(1), node->coord(2)));
  }

  pcl::toROSMsg(nodes_pcl, nodes_pcl_ros);
  nodes_pcl_ros.header.frame_id = "map";
  vis_node_pub.publish(nodes_pcl_ros);
}

void AStar::visConnections() {
  visualization_msgs::Marker connections;

  connections.header.frame_id = "map";
  connections.action = visualization_msgs::Marker::ADD;
  connections.pose.orientation.w = 1.0;
  connections.type = visualization_msgs::Marker::LINE_LIST;
  connections.scale.x = 0.1;
  connections.color.g = 1.0;
  connections.color.a = 0.5;

  int num_nodes = node_pool.size();
  for (int i = 0; i < num_nodes; i++) {
    NodePtr cur_node = node_pool.at(i);
    int num_connect_nodes = cur_node->connected_node.size();
    // if (cur_node->isGate) {
    //   ROS_WARN("(%f, %f, %f) is a GATE node having %d connections",
    //   cur_node->coord(0),
    //            cur_node->coord(1), cur_node->coord(2), num_connect_nodes);
    // } else {
    //   ROS_WARN("(%f, %f, %f) is a CENTER node having %d connections",
    //   cur_node->coord(0),
    //            cur_node->coord(1), cur_node->coord(2), num_connect_nodes);
    // }
    for (int j = 0; j < num_connect_nodes; j++) {
      NodePtr connect_node = cur_node->connected_node.at(j);

      geometry_msgs::Point p1, p2;
      p1.x = cur_node->coord(0);
      p1.y = cur_node->coord(1);
      p1.z = cur_node->coord(2);
      p2.x = connect_node->coord(0);
      p2.y = connect_node->coord(1);
      p2.z = connect_node->coord(2);

      connections.points.push_back(p1);
      connections.points.push_back(p2);
    }
  }

  vis_connections_pub.publish(connections);
}

void AStar::visParents() {
  visualization_msgs::Marker marker;

  marker.header.frame_id = "map";
  marker.action = visualization_msgs::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.type = visualization_msgs::Marker::LINE_LIST;
  marker.scale.x = 0.1;
  marker.color.g = 1.0;
  marker.color.a = 0.5;

  int num_nodes = node_pool.size();
  for (int i = 0; i < num_nodes; i++) {
    NodePtr cur_node = node_pool.at(i);
    NodePtr parent = cur_node->parent;

    if (parent == NULL) continue;

    geometry_msgs::Point p1, p2;
    p1.x = cur_node->coord(0);
    p1.y = cur_node->coord(1);
    p1.z = cur_node->coord(2);
    p2.x = parent->coord(0);
    p2.y = parent->coord(1);
    p2.z = parent->coord(2);

    marker.points.push_back(p1);
    marker.points.push_back(p2);
  }

  vis_parents_pub.publish(marker);
}

}  // namespace a_star
