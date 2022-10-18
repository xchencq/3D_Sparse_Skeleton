#include <math.h>
#include <Eigen/Eigen>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "../utils/quickhull/QuickHull.hpp"
#include "utils/A_star.h"
#include "utils/kdtree.h"

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/search/impl/kdtree.hpp>

#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/PositionCommand.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include "polygon_generation/backward.hpp"
#include "polygon_generation/data_type_3D.h"

inline bool compareFrontier(FrontierPtr f1, FrontierPtr f2) {
  return f1->facets.size() > f2->facets.size();
}

class SkeletonFinder {
 private:
  vector<shared_ptr<pcl::search::KdTree<pcl::PointXYZ>>> kdtreesForPolys;
  pcl::search::KdTree<pcl::PointXYZ> kdtreeForRawMap;
  pcl::search::KdTree<pcl::PointXYZ> kdtreeForVisMap;
  pcl::search::KdTree<pcl::PointXYZ> kdtreeForTestMap;
  pcl::search::KdTree<pcl::PointXYZ> kdtreeForNodes;

  // Map used in the main search process,
  // including obstacles, bounding box(optional), polyhedron facets
  pcl::PointCloud<pcl::PointXYZ> map_pcl;
  // Map used in closing loops,
  // including obstacles, bounding box
  pcl::PointCloud<pcl::PointXYZ> raw_map_pcl;
  // Map used in visualization,
  // including obstacles, excluding bounding box
  pcl::PointCloud<pcl::PointXYZ> vis_map_pcl;

  kdtree *kdTree_;  // dynamic light-weight Kd-tree, for organizing the nodes in
                    // the skeleton

  a_star::AStar path_finder;

  // deque<NodePtr> pending_nodes;
  deque<FrontierPtr> pending_frontiers;
  deque<FrontierPtr> loop_candidate_frontiers;
  vector<NodePtr> NodeList;
  vector<NodePtr> center_NodeList;
  vector<Eigen::Vector3d> sample_directions;
  vector<vector<Eigen::Vector3d>> bw_facets_directions;

  vector<int> pointIdxRadiusSearch;
  vector<float> pointRadiusSquaredDistance;
  vector<int> pointIdxRadiusSearchForRawMap;
  vector<float> pointRadiusSquaredDistanceForRawMap;
  vector<int> pointIdxRadiusSearchForVisMap;
  vector<float> pointRadiusSquaredDistanceForVisMap;
  vector<int> pointIdxRadiusSearchForNodes;
  vector<float> pointRadiusSquaredDistanceForNodes;

  Eigen::MatrixXd Path;
  Eigen::VectorXd Radius;
  vector<Eigen::Vector3d> path;

  ros::Subscriber map_sub;
  ros::Publisher vis_nodes_pub, vis_black_vertices_pub, vis_white_vertices_pub,
      vis_grey_vertices_pub, vis_map_pub, vis_polygons_pub, vis_frontiers_pub, vis_cur_frontier_pub,
      vis_connections_pub, vis_spheres_pub, vis_path_pub, vis_start_pub;
  pcl::PointCloud<pcl::PointXYZ> nodes_pcl, black_vertices_pcl, white_vertices_pcl,
      grey_vertices_pcl;
  sensor_msgs::PointCloud2 nodes_pcl_ros, black_vertices_pcl_ros, white_vertices_pcl_ros,
      grey_vertices_pcl_ros, map_pcl_ros, raw_map_pcl_ros, vis_map_pcl_ros;

  bool callback_executed = false;

  /* ------------------------ Param ----------------------- */
  // Map representation
  // 0: point cloud; 1: occupancy map
  int _map_representation;
  // Whether the map is in simulation
  bool _is_simulation;
  // An edge will be considered as a frontier if:
  // the dist to its nearest point exceeds this threshold
  double _frontier_creation_threshold;
  // Jump frontier
  double _frontier_jump_threshold;
  // Facets will be split into diff frontier if the angle between exceeds this threshold
  double _frontier_split_threshold;
  // A flowback will be created if number of contact vertices exceeds this threshold
  int _min_flowback_creation_threshold;
  // A flowback will not be created if the radius of contact vertices is below this threshold
  double _min_flowback_creation_radius_threshold;
  // A node will be discarded if its average vertex-center distance is below
  // this threshold
  double _min_node_radius;
  // A point on the ray will be considered as hit the pcl if:
  // the dist to its nearest point is below this margin
  // search_margin > sqrt((resolution/2)^2 + (raycast_step/2)^2)
  double _search_margin;
  // A ray will be discarded if length exceeds this max
  double _max_ray_length;
  // A new node will be set at the midpoint if length exceeds this max
  double _max_expansion_ray_length;
  // A node will be absorbed if difference of distance to floor with its parent exceeds this limit
  double _max_height_diff;
  // Number of sampings on the unit sphere
  int _sampling_density;
  // Max number of facets grouped in a frontier
  int _max_facets_grouped;
  // Resolution for map, raycast,
  double _resolution;
  // Visualization
  double _truncated_z_high;
  double _truncated_z_low;
  // Bounding box
  double _x_min, _x_max, _y_min, _y_max, _z_min, _z_max;
  double _start_x, _start_y, _start_z;
  double _path_start_x, _path_start_y, _path_start_z;
  double _path_target_x, _path_target_y, _path_target_z;

  /* ------------------ Development Tune ------------------ */
  bool _debug_mode;
  bool _bad_loop;

  // Visualize only the final result or the expansion process
  bool _visualize_final_result_only;
  // Visualize all or only the newest polyhedron
  bool _visualize_all;
  // Visualize outwards normal for each frontier
  bool _visualize_outwards_normal;
  // Visualize neighborhood facets for each frontier
  bool _visualize_nbhd_facets;
  // Visualize only_black polygon or black_and_white polygon
  bool _visualize_black_polygon;

  // Timers
  double genBlackAndWhiteVertices_timing = 0;
  double convex_hull_timing = 0;
  double centralizeNodePos_timing = 0;
  double findFlowBack_timing = 0;
  double identifyFacets_timing = 0;
  double identifyFrontiers_timing = 0;
  double addFacetsToPcl_timing = 0;
  double verifyFrontier_timing = 0;
  double processFrontier_timing = 0;

  /* -------------------- Sampling Use -------------------- */
  random_device rd;
  default_random_engine eng;

  // random distribution for generating samples inside a unit circle
  uniform_real_distribution<double> rand_rho = uniform_real_distribution<double>(0.0, 1.0);
  uniform_real_distribution<double> rand_phi = uniform_real_distribution<double>(0.0, 2 * M_PI);

  // basic random distributions for generating samples, in all feasible regions
  uniform_real_distribution<double> rand_x, rand_y, rand_z, rand_bias;
  // random distribution, especially for generating samples inside the local
  // map's boundary
  uniform_real_distribution<double> rand_x_in, rand_y_in, rand_z_in;

 public:
  SkeletonFinder();
  ~SkeletonFinder();

  /* set-up functions */
  void reset();
  void init(ros::NodeHandle &n);
  void setParam(ros::NodeHandle &n);
  void setStartPt(Eigen::Vector3d startPt);

  /* main function entries */
  void skeletonExpansion(Eigen::Vector3d startPt);
  bool initNode(NodePtr curNodePtr);
  void checkFacetClearance(NodePtr curNodePtr);
  bool checkPosConnected(NodePtr node_ptr, Eigen::Vector3d pos);
  void genBlackAndWhiteVertices(NodePtr nodePtr);
  void genSamplesOnUnitSphere();
  pair<Eigen::Vector3d, int> raycast(Eigen::Vector3d ray_source, Eigen::Vector3d direction,
                                     double cut_off_length);
  pair<Eigen::Vector3d, int> raycastOnRawMap(Eigen::Vector3d ray_source, Eigen::Vector3d direction,
                                             double cut_off_length);
  void centralizeNodePos(NodePtr node);
  void identifyBwFacets();
  void identifyFacets(NodePtr node);
  void findNbhdFacets(vector<FacetPtr> facets);
  void identifyFrontiers(NodePtr node);
  vector<FrontierPtr> splitFrontier(NodePtr node, vector<FacetPtr> group_facets);
  vector<FacetPtr> findGroupFacetsFromVertices(NodePtr node, vector<VertexPtr> group_bv);
  void findFlowBack(NodePtr node);
  void combineFrontier(NodePtr node);
  bool initFrontier(FrontierPtr frontier);
  bool checkFacetsCombine(FacetPtr f1, FacetPtr f2);
  bool checkPtInPolyhedron(NodePtr node, Eigen::Vector3d pt);
  void verifyFrontier(FrontierPtr ftr_ptr);
  void addFacetsToPcl(NodePtr nodePtr);
  bool processFrontier(FrontierPtr curFtrPtr);
  vector<Eigen::Vector3d> findPath(Eigen::Vector3d start, Eigen::Vector3d target);
  vector<Eigen::Vector3d> findPathByAStar(Eigen::Vector3d start, Eigen::Vector3d target);

  bool isSamePos(Eigen::Vector3d pos1, Eigen::Vector3d pos2);
  bool checkFloor(NodePtr node);
  int onCeilOrFloor(Eigen::Vector3d p);
  bool facetOnCeilOrFloor(FacetPtr f);
  bool checkSegClear(Eigen::Vector3d pos1, Eigen::Vector3d pos2);
  double getNodeRadius(NodePtr curNodePtr);
  double getVerticesRadius(vector<Eigen::Vector3d> vertices);
  VertexPtr getVertexFromDire(NodePtr node, Eigen::Vector3d dire);
  FacetPtr getCommonInvalidNbhdFacet(FrontierPtr f1, FrontierPtr f2, Eigen::Vector3d pt);
  vector<Eigen::Vector3d> getCommonVertices(FrontierPtr f1, FrontierPtr f2);
  bool checkPtOnFrontier(FrontierPtr ftr_ptr, Eigen::Vector3d pt);
  bool checkPtOnFacet(FacetPtr facet, Eigen::Vector3d pt);
  FrontierPtr findFlowBackFrontier(Eigen::Vector3d hit_on_pcl, int index);
  FacetPtr findFlowBackFacet(Eigen::Vector3d pos, int index);
  bool checkWithinBbx(Eigen::Vector3d pos);
  void addBbxToMap(pcl::PointCloud<pcl::PointXYZ> &map);

  void mapCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);

  /* operations on the tree */
  void recordNode(NodePtr new_node);
  void treeDestruct();

  /* utility functions */
  inline double getDis(const NodePtr node1, const NodePtr node2);
  inline double getDis(const NodePtr node1, const Eigen::Vector3d &pt);
  inline double getDis(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);
  inline double getDis(const Eigen::Vector3i &p1, const Eigen::Vector3i &p2);
  // inline Eigen::Vector3d genSample();
  inline pair<double, int> radiusSearch(Eigen::Vector3d &pt);
  inline double radiusSearchOnRawMap(Eigen::Vector3d &pt);
  inline double radiusSearchOnVisMap(Eigen::Vector3d &pt);
  bool collisionCheck(Eigen::Vector3d search_pt, double threshold);

  // inline NodePtr genNewNode(Eigen::Vector3d &pt);
  // NodePtr pending_nodes_pop_front() {
  //   NodePtr curNodePtr = pending_nodes.front();
  //   pending_nodes.pop_front();
  //   return curNodePtr;
  // };
  FrontierPtr pendingFrontiersPopFront() {
    FrontierPtr curFrontierPtr = pending_frontiers.front();
    pending_frontiers.pop_front();
    return curFrontierPtr;
  };

  /* data return */
  pair<Eigen::MatrixXd, Eigen::VectorXd> getPath() { return make_pair(Path, Radius); };

  /* -------------------- visualization ------------------- */
  void visualization();
  void visNodesAndVertices();
  void visMap();
  void visPolygons();
  void visFrontiers();
  void visCurrentFrontier(FrontierPtr ftr);
  void visConnections();
  void visSphere(Eigen::Vector3d pos, double radius);
  void visPath();
  void visStart();
};