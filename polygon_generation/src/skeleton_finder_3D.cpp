#include "polygon_generation/skeleton_finder_3D.h"
#include <algorithm>

using namespace Eigen;
using namespace std;

namespace backward {
backward::SignalHandling sh;
}

SkeletonFinder::SkeletonFinder() {}

SkeletonFinder::~SkeletonFinder() {}

void SkeletonFinder::setParam(ros::NodeHandle &n) {
  n.param("mapBoundary/x_min", _x_min, -10.0);
  n.param("mapBoundary/x_max", _x_max, 10.0);
  n.param("mapBoundary/y_min", _y_min, -10.0);
  n.param("mapBoundary/y_max", _y_max, 10.0);
  n.param("mapBoundary/z_min", _z_min, 0.0);
  n.param("mapBoundary/z_max", _z_max, 3.0);
  n.param("start_x", _start_x, 0.0);
  n.param("start_y", _start_y, 0.0);
  n.param("start_z", _start_z, 0.0);
  n.param("path_start_x", _path_start_x, 0.0);
  n.param("path_start_y", _path_start_y, 0.0);
  n.param("path_start_z", _path_start_z, 0.0);
  n.param("path_target_x", _path_target_x, 0.0);
  n.param("path_target_y", _path_target_y, 0.0);
  n.param("path_target_z", _path_target_z, 0.0);

  n.param("map_representation", _map_representation, 0);
  n.param("is_simulation", _is_simulation, true);
  n.param("frontier_creation_threshold", _frontier_creation_threshold, 2.0);
  n.param("frontier_jump_threshold", _frontier_jump_threshold, 2.0);
  n.param("frontier_split_threshold", _frontier_split_threshold, 2.0);
  n.param("min_loop_creation_threshold", _min_flowback_creation_threshold, 5);
  n.param("min_flowback_creation_radius_threshold",
          _min_flowback_creation_radius_threshold, 1.0);

  n.param("min_node_radius", _min_node_radius, 2.0);
  n.param("search_margin", _search_margin, 0.2);
  n.param("max_ray_length", _max_ray_length, 8.0);
  n.param("max_expansion_ray_length", _max_expansion_ray_length, 5.0);
  n.param("max_height_diff", _max_height_diff, 5.0);
  n.param("sampling_density", _sampling_density, 16);
  n.param("max_facets_grouped", _max_facets_grouped, 6);
  n.param("resolution", _resolution, 0.2);
  n.param("truncated_z_high", _truncated_z_high, 2.5);
  n.param("truncated_z_low", _truncated_z_low, 0.0);

  n.param("debug_mode", _debug_mode, true);
  n.param("bad_loop", _bad_loop, true);

  n.param("visualize_final_result_only", _visualize_final_result_only, true);
  n.param("visualize_all", _visualize_all, true);
  n.param("visualize_outwards_normal", _visualize_outwards_normal, false);
  n.param("visualize_nbhd_frontiers", _visualize_nbhd_facets, false);
  n.param("visualize_black_polygon", _visualize_black_polygon, false);
}

void SkeletonFinder::reset() {
  treeDestruct();
  NodeList.clear();
}

void SkeletonFinder::setStartPt(Vector3d startPt) {
  NodePtr start = new Node(startPt, NULL);
  initNode(start);
}

inline double SkeletonFinder::getDis(const NodePtr node1, const NodePtr node2) {
  return sqrt(pow(node1->coord(0) - node2->coord(0), 2) +
              pow(node1->coord(1) - node2->coord(1), 2) +
              pow(node1->coord(2) - node2->coord(2), 2));
}

inline double SkeletonFinder::getDis(const NodePtr node1, const Vector3d &pt) {
  return sqrt(pow(node1->coord(0) - pt(0), 2) +
              pow(node1->coord(1) - pt(1), 2) +
              pow(node1->coord(2) - pt(2), 2));
}

inline double SkeletonFinder::getDis(const Vector3d &p1, const Vector3d &p2) {
  return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2) +
              pow(p1(2) - p2(2), 2));
}

inline double SkeletonFinder::getDis(const Vector3i &p1, const Vector3i &p2) {
  return sqrt((double)(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2) +
                       pow(p1(2) - p2(2), 2)));
}

inline pair<double, int> SkeletonFinder::radiusSearch(Vector3d &search_Pt) {
  double min_dis = radiusSearchOnRawMap(search_Pt);
  int min_dis_node_index = -1;

  if (min_dis < _search_margin) {
    pair<double, int> return_pair(min_dis, min_dis_node_index);
    return return_pair;
  }

  // if (collisionCheck(search_Pt, _search_margin)) {
  //   pair<double, int> return_pair(min_dis, min_dis_node_index);
  //   return return_pair;
  // }

  pcl::PointXYZ searchPoint;
  searchPoint.x = search_Pt(0);
  searchPoint.y = search_Pt(1);
  searchPoint.z = search_Pt(2);

  for (NodePtr node : NodeList) {
    if (node->rollbacked || node->isGate)
      continue;

    if (getDis(node->original_coord, search_Pt) > _max_ray_length + min_dis)
      continue;

    pointIdxRadiusSearch.clear();
    pointRadiusSquaredDistance.clear();

    kdtreesForPolys.at(node->index)
        ->nearestKSearch(searchPoint, 1, pointIdxRadiusSearchForRawMap,
                         pointRadiusSquaredDistanceForRawMap);
    double radius = sqrt(pointRadiusSquaredDistanceForRawMap[0]);

    if (radius < min_dis) {
      min_dis = radius;
      min_dis_node_index = node->index;
    }
  }

  pair<double, int> return_pair(min_dis, min_dis_node_index);
  return return_pair;
}

inline double SkeletonFinder::radiusSearchOnRawMap(Vector3d &search_Pt) {
  //   if (getDis(search_Pt, start_pt) > sample_range + max_radius)
  //     return max_radius - search_margin;

  pcl::PointXYZ searchPoint;
  searchPoint.x = search_Pt(0);
  searchPoint.y = search_Pt(1);
  searchPoint.z = search_Pt(2);

  pointIdxRadiusSearchForRawMap.clear();
  pointRadiusSquaredDistanceForRawMap.clear();

  kdtreeForRawMap.nearestKSearch(searchPoint, 1, pointIdxRadiusSearchForRawMap,
                                 pointRadiusSquaredDistanceForRawMap);
  double radius = sqrt(pointRadiusSquaredDistanceForRawMap[0]);
  //   return min(radius, double(max_radius));
  return radius;
}

bool SkeletonFinder::collisionCheck(Vector3d search_pt, double threshold) {
  // Point cloud map
  if (_map_representation == 0) {
    double dis = radiusSearchOnRawMap(search_pt);
    if (dis < threshold)
      return true;
    else
      return false;
  }
  // Won't reach
  return true;
}
void SkeletonFinder::recordNode(NodePtr new_node) {
  NodeList.push_back(new_node);
  nodes_pcl.points.push_back(pcl::PointXYZ(
      new_node->coord(0), new_node->coord(1), new_node->coord(2)));

  if (!new_node->isGate) {
    int index = center_NodeList.size();
    new_node->index = index;
    center_NodeList.push_back(new_node);
  }

  // ROS_INFO("New node recorded.");
  // ROS_INFO("Current node list size: %i", NodeList.size());
}

void SkeletonFinder::treeDestruct() {
  kd_free(kdTree_);

  for (int i = 0; i < int(NodeList.size()); i++) {
    NodePtr ptr = NodeList[i];
    delete ptr;
  }
}

/* ------------------------------------------------------ */

void SkeletonFinder::init(ros::NodeHandle &n) {
  /* --------------------- Set params --------------------- */
  setParam(n);
  path_finder.initHandler(n);

  map_sub = n.subscribe("/map/map", 1, &SkeletonFinder::mapCallBack, this);

  vis_nodes_pub =
      n.advertise<sensor_msgs::PointCloud2>("/skeleton_finder_3D/nodes", 1);

  vis_black_vertices_pub = n.advertise<sensor_msgs::PointCloud2>(
      "/skeleton_finder_3D/black_vertices", 1);

  vis_white_vertices_pub = n.advertise<sensor_msgs::PointCloud2>(
      "/skeleton_finder_3D/white_vertices", 1);

  vis_grey_vertices_pub = n.advertise<sensor_msgs::PointCloud2>(
      "/skeleton_finder_3D/grey_vertices", 1);

  vis_map_pub =
      n.advertise<sensor_msgs::PointCloud2>("/skeleton_finder_3D/map", 1);

  vis_polygons_pub = n.advertise<visualization_msgs::Marker>(
      "/skeleton_finder_3D/polygons", 1);

  vis_frontiers_pub =
      n.advertise<visualization_msgs::Marker>("/skeleton_finder_3D/facets", 1);

  vis_cur_frontier_pub = n.advertise<visualization_msgs::Marker>(
      "/skeleton_finder_3D/frontier", 1);

  vis_connections_pub = n.advertise<visualization_msgs::Marker>(
      "/skeleton_finder_3D/connections", 1);

  vis_spheres_pub =
      n.advertise<visualization_msgs::Marker>("/skeleton_finder_3D/spheres", 1);

  vis_path_pub =
      n.advertise<visualization_msgs::Marker>("/skeleton_finder_3D/path", 1);

  vis_start_pub =
      n.advertise<visualization_msgs::Marker>("/skeleton_finder_3D/start", 1);
}

void SkeletonFinder::mapCallBack(
    const sensor_msgs::PointCloud2 &pointcloud_map) {
  if (callback_executed)
    return;

  callback_executed = true;

  pcl::fromROSMsg(pointcloud_map, map_pcl);
  pcl::fromROSMsg(pointcloud_map, raw_map_pcl);

  if (map_pcl.points.size() == 0) {
    ROS_WARN("Map pointcloud is empty!");
    return;
  }
  ROS_INFO("map_pcl size: %d", (int)map_pcl.points.size());
  for (int i = 0; i < (int)map_pcl.points.size(); i++) {
    if (map_pcl.points[i].z < 0.2)
      continue;
    // if (map_pcl.points[i].z < _truncated_z_low || map_pcl.points[i].z >
    // _truncated_z_high) continue;
    vis_map_pcl.points.push_back(map_pcl.points[i]);
  }
  vis_map_pcl.width = vis_map_pcl.points.size();
  vis_map_pcl.height = 1;
  // ROS_INFO("vis_map size: %d", vis_map_pcl.points.size());

  addBbxToMap(raw_map_pcl);
  addBbxToMap(map_pcl);

  // Point cloud map
  if (_map_representation == 0) {
    kdtreeForRawMap.setInputCloud(raw_map_pcl.makeShared());
  }

  visualization();

  ROS_ERROR("Wait key to start...");
  getchar();

  Eigen::Vector3d start;
  start << _start_x, _start_y, _start_z;

  ROS_INFO("Generating skeleton...");
  ros::Time begin = ros::Time::now();
  skeletonExpansion(start);
  ros::Time finish = ros::Time::now();

  ROS_INFO("Expansion finished.");
  double timing = (finish - begin).toSec();

  ROS_INFO("Expansion time: %f", timing);
  // ROS_ERROR("genBlackAndWhiteVertices_timing takes up: %f (%f percents)",
  //           genBlackAndWhiteVertices_timing, genBlackAndWhiteVertices_timing
  //           / timing);
  // ROS_ERROR("convex_hull_timing takes up: %f (%f percents)",
  // convex_hull_timing,
  //           convex_hull_timing / timing);
  // ROS_ERROR("centralizeNodePos_timing takes up: %f (%f percents)",
  // centralizeNodePos_timing,
  //           centralizeNodePos_timing / timing);
  // ROS_ERROR("findFlowBack_timing takes up: %f (%f percents)",
  // findFlowBack_timing,
  //           findFlowBack_timing / timing);
  // ROS_ERROR("identifyFacets_timing takes up: %f (%f percents)",
  // identifyFacets_timing,
  //           identifyFacets_timing / timing);
  // ROS_ERROR("identifyFrontiers_timing takes up: %f (%f percents)",
  // identifyFrontiers_timing,
  //           identifyFrontiers_timing / timing);
  // ROS_ERROR("addFacetsToPcl_timing takes up: %f (%f percents)",
  // addFacetsToPcl_timing,
  //           addFacetsToPcl_timing / timing);
  // ROS_ERROR("verifyFrontier_timing takes up: %f (%f percents)",
  // verifyFrontier_timing,
  //           verifyFrontier_timing / timing);
  // ROS_ERROR("processFrontier_timing takes up: %f (%f percents)",
  // processFrontier_timing,
  //           processFrontier_timing / timing);

  visualization();

  ROS_ERROR("Wait key to start finding path...");
  getchar();

  // Eigen::Vector3d start_query(-28, -28, 1);
  // Eigen::Vector3d target_query(28, 28, 1);
  Eigen::Vector3d start_query(_path_start_x, _path_start_y, _path_start_z);
  Eigen::Vector3d target_query(_path_target_x, _path_target_y, _path_target_z);
  begin = ros::Time::now();
  path = findPath(start_query, target_query);
  finish = ros::Time::now();
  if (path.empty()) {
    ROS_ERROR("Find path failed!");
  } else {
    visPath();
    ROS_WARN("Path finding time: %f", (finish - begin).toSec() * 1000);
    double path_length = path_finder.pathLength(path);
    ROS_WARN("Path length: %f", path_length);
  }
}

void SkeletonFinder::skeletonExpansion(Eigen::Vector3d startPt) {
  genSamplesOnUnitSphere();
  identifyBwFacets();
  setStartPt(startPt);
  // ROS_INFO("Start point set!");

  FrontierPtr cur_frontier;
  while (!pending_frontiers.empty()) {
    ros::Time begin = ros::Time::now();
    // visualization();
    // cout << "Wait for key to expand next frontier..." << endl;
    // getchar();

    cur_frontier = pendingFrontiersPopFront();
    if (cur_frontier->deleted)
      continue;
    if (cur_frontier == NULL)
      continue;

    // ROS_ERROR("Poped a new frontier: (%f, %f, %f)",
    // cur_frontier->proj_center(0),
    //           cur_frontier->proj_center(1), cur_frontier->proj_center(2));

    verifyFrontier(cur_frontier);
    // ROS_INFO("verifyFrontier trail 1: %d", cur_frontier->valid);
    // visCurrentFrontier(cur_frontier);
    // getchar();

    if (!cur_frontier->valid) {
      // ROS_INFO("Try... Change to use proj_facet_center");
      Eigen::Vector3d prev_proj_center = cur_frontier->proj_center;
      Eigen::Vector3d prev_normal = cur_frontier->outwards_unit_normal;
      cur_frontier->proj_center = cur_frontier->proj_facet->center;
      verifyFrontier(cur_frontier);
      // ROS_INFO("verifyFrontier trail 2: %d", cur_frontier->valid);
      // visCurrentFrontier(cur_frontier);
      // getchar();

      if (!cur_frontier->valid) {
        // ROS_INFO("Try... Change to use proj_facet_normal");
        cur_frontier->proj_center = prev_proj_center;
        cur_frontier->outwards_unit_normal =
            cur_frontier->proj_facet->outwards_unit_normal;
        verifyFrontier(cur_frontier);
        // ROS_INFO("verifyFrontier trail 3: %d", cur_frontier->valid);
        // visCurrentFrontier(cur_frontier);
        // getchar();

        if (!cur_frontier->valid) {
          // ROS_INFO("Try... Change to use proj_facet_center &&
          // proj_facet_normal");

          cur_frontier->proj_center = cur_frontier->proj_facet->center;
          // cur_frontier->outwards_unit_normal =
          // cur_frontier->proj_facet_normal; // already set
          verifyFrontier(cur_frontier);
          // ROS_INFO("verifyFrontier trail 4: %d", cur_frontier->valid);
          // visCurrentFrontier(cur_frontier);
          // getchar();

          if (!cur_frontier->valid) {
            cur_frontier->outwards_unit_normal = prev_normal;
            for (FacetPtr candidate_facet :
                 cur_frontier->proj_facet->nbhd_facets) {
              // ROS_INFO("Try... Change to use nbhd facet_center");
              cur_frontier->proj_center = candidate_facet->center;
              verifyFrontier(cur_frontier);
              // ROS_INFO("verifyFrontier trail 5: %d", cur_frontier->valid);
              // visCurrentFrontier(cur_frontier);
              // getchar();

              if (!cur_frontier->valid) {
                cur_frontier->outwards_unit_normal =
                    candidate_facet->outwards_unit_normal;
                verifyFrontier(cur_frontier);
                // ROS_INFO("verifyFrontier trail 6: %d", cur_frontier->valid);
                // visCurrentFrontier(cur_frontier);
                // getchar();
                if (cur_frontier->valid)
                  break;
              } else {
                break;
              }
            }
          }
        }
      }
    }
    ros::Time finish = ros::Time::now();
    verifyFrontier_timing += (finish - begin).toSec();
    begin = ros::Time::now();
    if (cur_frontier->valid) {
      // ROS_INFO_STREAM("Frontier next_node_pos: " <<
      // cur_frontier->next_node_pos.transpose());
      if (!processFrontier(cur_frontier)) {
        // ROS_INFO("processFrontier fails.");
      } else {
        if (!_visualize_final_result_only)
          visualization();
        if (_debug_mode) {
          ROS_WARN("Expanded a valid frontier.");
          visualization();
          getchar();
        }
      }
    } else {
      if (_debug_mode) {
        ROS_WARN("Frontier not valid!");
        // cout << "Frontier center:" << cur_frontier->proj_center.transpose()
        // << endl; visualization(); getchar();
      }
    }
    finish = ros::Time::now();
    processFrontier_timing += (finish - begin).toSec();
  }
}

bool SkeletonFinder::initNode(NodePtr curNodePtr) {
  // ROS_ERROR("Start init node");
  ros::Time begin = ros::Time::now();
  ros::Time prev = begin;
  ros::Time curr;

  if (!checkWithinBbx(curNodePtr->coord)) {
    return false;
  }
  if (!checkFloor(curNodePtr)) {
    // ROS_ERROR("No floor! Absort.");
    return false;
  }

  // if(curNodePtr->coord(1) < -18) return false;
  if (!curNodePtr->isGate) {
    // ROS_ERROR("Initing node pos:(%f, %f, %f)", curNodePtr->coord(0),
    // curNodePtr->coord(1),
    //           curNodePtr->coord(2));

    genBlackAndWhiteVertices(curNodePtr);

    curr = ros::Time::now();
    genBlackAndWhiteVertices_timing += (curr - prev).toSec();
    prev = curr;

    if (curNodePtr->black_vertices.size() < 4) {
      // ROS_WARN("Found only %d vertices for the current node!",
      // curNodePtr->black_vertices.size());
      return false;
    }

    centralizeNodePos(curNodePtr);

    curr = ros::Time::now();
    centralizeNodePos_timing += (curr - prev).toSec();
    prev = curr;

    double radius = getNodeRadius(curNodePtr);
    if (radius < _min_node_radius && curNodePtr->white_vertices.empty()) {
      // ROS_WARN("This node radius is too small!");
      // ROS_WARN("Node radius: %f", radius);
      return false;
    }

    // Absorb node inside node
    if (curNodePtr->seed_frontier != NULL) {
      bool diff_ind = false;
      int ind = curNodePtr->black_vertices.at(0)->collision_node_index;
      for (VertexPtr v : curNodePtr->black_vertices) {
        if (v->collision_node_index != ind) {
          diff_ind = true;
          break;
        }
      }
      if (!diff_ind)
        return false;
    }

    findFlowBack(curNodePtr);

    curr = ros::Time::now();
    findFlowBack_timing += (curr - prev).toSec();
    prev = curr;

    identifyFacets(curNodePtr);

    curr = ros::Time::now();
    identifyFacets_timing += (curr - prev).toSec();
    prev = curr;

    identifyFrontiers(curNodePtr);

    curr = ros::Time::now();
    identifyFrontiers_timing += (curr - prev).toSec();
    prev = curr;

    addFacetsToPcl(curNodePtr);

    curr = ros::Time::now();
    addFacetsToPcl_timing += (curr - prev).toSec();

    // visualization();
  }

  // curNodePtr->clearance = radiusSearchOnRawMap(curNodePtr->coord);
  // visSphere(curNodePtr->coord, curNodePtr->clearance);
  recordNode(curNodePtr);
  return true;
}

void SkeletonFinder::genBlackAndWhiteVertices(NodePtr nodePtr) {
  // vector<VertexPtr> black;
  vector<Eigen::Vector3d>::iterator it;
  for (it = sample_directions.begin(); it != sample_directions.end(); it++) {
    int index = nodePtr->sampling_directions.size();
    nodePtr->sampling_directions.push_back(vec3((*it)(0), (*it)(1), (*it)(2)));

    pair<Vector3d, int> raycast_result =
        raycast(nodePtr->coord, *it, _max_ray_length);
    Eigen::Vector3d newVertex = raycast_result.first;
    if (raycast_result.second == -2) {
      newVertex += (*it) * _max_ray_length;
      VertexPtr new_white_vertex = new Vertex(newVertex, (*it), WHITE);
      new_white_vertex->sampling_dire_index = index;
      nodePtr->white_vertices.push_back(new_white_vertex);
    } else {
      VertexPtr new_black_vertex;
      new_black_vertex = new Vertex(newVertex, (*it), BLACK);
      new_black_vertex->collision_node_index = raycast_result.second;
      /*
      // if (radiusSearchOnRawMap(newVertex) > _search_margin) {
      if (raycast_result.second > -1) {
        new_black_vertex = new Vertex(newVertex, (*it), GREY);
        new_black_vertex->collision_node_index = raycast_result.second;
      } else {
        new_black_vertex = new Vertex(newVertex, (*it), BLACK);
      }
      */
      new_black_vertex->sampling_dire_index = index;
      new_black_vertex->dis_to_center =
          getDis(new_black_vertex->coord, nodePtr->coord);
      // black.push_back(new_black_vertex);
      nodePtr->black_vertices.push_back(new_black_vertex);
      // nodePtr->valid_sampling_directions.push_back(vec3((*it)(0), (*it)(1),
      // (*it)(2)));
    }
  }
  // double avg_black_length = 0;
  // int black_count = 0;
  // for(VertexPtr v : black){
  //   avg_black_length = (avg_black_length * black_count + v->dis_to_center) /
  //   (black_count + 1); black_count++;
  // }
  // for(VertexPtr v : black){
  //   if(v->dis_to_center > 2 * avg_black_length){
  //     v->type = WHITE;
  //     nodePtr->white_vertices.push_back(v);
  //   }
  //   else{
  //     nodePtr->black_vertices.push_back(v);
  //   }
  // }
}

void SkeletonFinder::genSamplesOnUnitSphere() {
  // Fibonicci sphere
  double phi = M_PI * (3 - sqrt(5));
  double x, y, z, radius, theta;

  for (int i = 0; i < _sampling_density; i++) {
    y = 1 - 2 * ((float)i / (float)(_sampling_density - 1));
    radius = sqrt(1 - y * y);
    theta = phi * i;
    x = cos(theta) * radius;
    z = sin(theta) * radius;

    Eigen::Vector3d sample;
    sample << x, y, z;
    sample_directions.push_back(sample);
  }
}

// -2: collision not found within cut_off_length
// -1: collision is with map
// >=0: collision is with the node of that index
pair<Vector3d, int> SkeletonFinder::raycast(Vector3d ray_source,
                                            Vector3d direction,
                                            double cut_off_length) {
  double clearance = radiusSearch(ray_source).first;
  if (clearance > cut_off_length) {
    pair<Vector3d, int> return_pair(ray_source, -2);
    return return_pair;
  } else {
    // Eigen::Vector3d offset = _resolution * direction;
    Eigen::Vector3d current_pos = ray_source + clearance * direction;
    // int cnt = 0;
    double length = clearance;
    // ROS_ERROR("Ray_source: %f, %f, %f", ray_source(0), ray_source(1),
    // ray_source(2)); ROS_ERROR("Direction:%f, %f, %f", direction(0),
    // direction(1), direction(2));

    // while (cnt * _resolution <= cut_off_length) {
    // while (cnt * _resolution <= cut_off_length - clearance) {
    while (length <= cut_off_length) {
      pair<double, int> rs = radiusSearch(current_pos);
      double radius = rs.first;
      // ROS_WARN("radiusSearch value: %f", radius);

      if (radius < _search_margin) {
        pair<Vector3d, int> return_pair(current_pos, rs.second);
        return return_pair;
      }

      // current_pos += offset;
      // cnt++;
      current_pos += radius * direction;
      length += radius;
    }

    pair<Vector3d, int> return_pair(ray_source, -2);
    return return_pair;
  }
}

// -2: collision not found within cut_off_length
// -1: collision is with map
pair<Vector3d, int> SkeletonFinder::raycastOnRawMap(Vector3d ray_source,
                                                    Vector3d direction,
                                                    double cut_off_length) {
  // Point cloud map
  if (_map_representation == 0) {
    double clearance = radiusSearchOnRawMap(ray_source);
    if (clearance > cut_off_length) {
      pair<Vector3d, int> return_pair(ray_source, -2);
      return return_pair;
    } else {
      Eigen::Vector3d current_pos = ray_source + clearance * direction;
      double length = clearance;

      while (length <= cut_off_length) {
        double radius = radiusSearchOnRawMap(current_pos);

        if (radius < _search_margin) {
          pair<Vector3d, int> return_pair(current_pos, -1);
          return return_pair;
        }
        current_pos += radius * direction;
        length += radius;
      }

      pair<Vector3d, int> return_pair(ray_source, -2);
      return return_pair;
    }
  }
  // Won't reach
  pair<Vector3d, int> temp(Vector3d::Zero(), 0);
  return temp;
}

void SkeletonFinder::centralizeNodePos(NodePtr node) {
  int cnt = 0;
  Vector3d sum = Vector3d::Zero();

  vector<VertexPtr>::iterator it;
  for (it = node->black_vertices.begin(); it != node->black_vertices.end();
       it++) {
    cnt++;
    sum = sum + (*it)->coord;
  }
  node->coord = sum / cnt;
}

void SkeletonFinder::identifyBwFacets() {
  // Mesh2 is for black and white polygon
  quickhull::QuickHull<double> qh;
  quickhull::HalfEdgeMesh<double, size_t> mesh2 = qh.getConvexHullAsMesh(
      &sample_directions[0](0), sample_directions.size(), true);

  for (auto &face : mesh2.m_faces) {
    quickhull::HalfEdgeMesh<double, size_t>::HalfEdge &halfedge =
        mesh2.m_halfEdges[face.m_halfEdgeIndex];

    vec3 vertex_quickhull;
    vector<Eigen::Vector3d> vertices_eigen;
    Eigen::Vector3d vertex_eigen;
    for (int i = 0; i < 3; i++) {
      vertex_quickhull = mesh2.m_vertices[halfedge.m_endVertex];
      vertex_eigen = Eigen::Vector3d::Zero();
      vertex_eigen << vertex_quickhull.x, vertex_quickhull.y,
          vertex_quickhull.z;
      vertices_eigen.push_back(vertex_eigen);
      halfedge = mesh2.m_halfEdges[halfedge.m_next];
    }

    bw_facets_directions.push_back(vertices_eigen);
  }
}

void SkeletonFinder::identifyFacets(NodePtr node) {
  for (vector<Eigen::Vector3d> facet_vertices : bw_facets_directions) {
    VertexPtr v1, v2, v3;
    v1 = getVertexFromDire(node, facet_vertices.at(0));
    v2 = getVertexFromDire(node, facet_vertices.at(1));
    v3 = getVertexFromDire(node, facet_vertices.at(2));
    v1->connected_vertices.push_back(v2);
    v2->connected_vertices.push_back(v3);
    v3->connected_vertices.push_back(v1);
  }
}

void SkeletonFinder::identifyFrontiers(NodePtr node) {
  vector<vector<VertexPtr>> bv_groups;
  int num_wv = node->white_vertices.size();
  for (int i = 0; i < num_wv; i++) {
    VertexPtr seed_wv = node->white_vertices.at(i);
    if (seed_wv->visited)
      continue;

    seed_wv->visited = true;
    vector<VertexPtr> group_bv;
    deque<VertexPtr> pending_wv;

    pending_wv.push_back(seed_wv);

    // ROS_ERROR("A new group of black points:");
    // while (!pending_wv.empty() && group_bv.size() < _max_facets_grouped) {
    while (!pending_wv.empty()) {
      VertexPtr v = pending_wv.front();
      v->visited = true;
      pending_wv.pop_front();

      for (VertexPtr v_nbhd : v->connected_vertices) {
        if (v_nbhd->type == WHITE) {
          if (v_nbhd->visited)
            continue;
          Eigen::Vector3d midpt = (v->coord + v_nbhd->coord) / 2;
          if (radiusSearch(midpt).first > 2 * _search_margin) {
            pending_wv.push_back(v_nbhd);
          }
        } else {
          // v_nbhd->visited = true;
          v_nbhd->critical = true;
          group_bv.push_back(v_nbhd);
          // ROS_WARN("(%f, %f, %f)", v_nbhd->coord(0), v_nbhd->coord(1),
          // v_nbhd->coord(2));
        }
      }
    }

    if (group_bv.size() < 3)
      continue;
    for (VertexPtr v : group_bv)
      v->visited = true;
    bv_groups.push_back(group_bv);
  }

  // Filter black vertices
  int num_groups = bv_groups.size();
  for (int i = 0; i < num_groups; i++) {
    double mean_length = 0;
    double tolerance;
    for (VertexPtr v : bv_groups.at(i)) {
      mean_length += v->dis_to_center;
    }
    mean_length /= bv_groups.at(i).size();
    tolerance = mean_length * 0.3;
    int longest_index = -1;
    int shortest_index = -1;
    double longest = 0;
    double shortest = 9999;
    int num_ver = bv_groups.at(i).size();
    for (int j = 0; j < num_ver; j++) {
      if (onCeilOrFloor(bv_groups.at(i).at(j)->coord) == 0)
        continue;
      double dis = bv_groups.at(i).at(j)->dis_to_center;
      if (dis - mean_length > tolerance) {
        if (dis > longest) {
          longest_index = j;
          longest = dis;
        }
      } else if (mean_length - dis > tolerance) {
        if (dis < shortest) {
          shortest_index = j;
          shortest = dis;
        }
      }
    }
    if (longest_index != -1) {
      bv_groups.at(i).at(longest_index)->type = GREY;
      bv_groups.at(i).at(longest_index)->critical = false;
    }
    if (shortest_index != -1) {
      bv_groups.at(i).at(shortest_index)->type = GREY;
      bv_groups.at(i).at(shortest_index)->critical = false;
    }
  }

  // Inflate critical black vertices
  for (int i = 0; i < num_groups; i++) {
    for (VertexPtr v : bv_groups.at(i)) {
      for (VertexPtr v_nbhd : v->connected_vertices) {
        if (v_nbhd->type == BLACK)
          v_nbhd->critical = true;
      }
    }
  }

  // Mesh1 is for only black polygon
  vector<vec3> bv_for_mesh;
  for (VertexPtr bv : node->black_vertices) {
    if (onCeilOrFloor(bv->coord) != 0 && !bv->critical)
      continue;
    bv_for_mesh.push_back(
        node->sampling_directions.at(bv->sampling_dire_index));
    // node->poly_vertices.push_back(bv);
    bv->critical = true;
  }

  quickhull::QuickHull<double> qh;
  quickhull::HalfEdgeMesh<double, size_t> mesh1 =
      qh.getConvexHullAsMesh(&(bv_for_mesh)[0].x, (bv_for_mesh).size(), true);

  for (auto &face : mesh1.m_faces) {
    vector<VertexPtr> vertices;

    quickhull::HalfEdgeMesh<double, size_t>::HalfEdge &halfedge =
        mesh1.m_halfEdges[face.m_halfEdgeIndex];

    for (int i = 0; i < 3; i++) {
      vec3 vertex_quickhull = mesh1.m_vertices[halfedge.m_endVertex];
      Eigen::Vector3d vertex_eigen;
      vertex_eigen << vertex_quickhull.x, vertex_quickhull.y,
          vertex_quickhull.z;
      vertices.push_back(getVertexFromDire(node, vertex_eigen));

      halfedge = mesh1.m_halfEdges[halfedge.m_next];
    }

    FacetPtr new_facet = new Facet(vertices, node);
    int ind = node->facets.size();
    node->facets.push_back(new_facet);
    new_facet->index = ind;
  }
  // ROS_ERROR("identifyFacets: %d", node->facets.size());

  // Calculate outwards normal for each facet
  int num_facet = node->facets.size();
  for (int i = 0; i < num_facet; i++) {
    FacetPtr facet_ptr = node->facets.at(i);

    Eigen::Vector3d v1 =
        facet_ptr->vertices.at(1)->coord - facet_ptr->vertices.at(0)->coord;
    Eigen::Vector3d v2 =
        facet_ptr->vertices.at(2)->coord - facet_ptr->vertices.at(0)->coord;
    Eigen::Vector3d candidate_normal = v1.cross(v2);
    candidate_normal.normalize();

    Eigen::Vector3d pt_to_judge = facet_ptr->center + candidate_normal;
    if (checkPtInPolyhedron(node, pt_to_judge))
      facet_ptr->outwards_unit_normal = -candidate_normal;
    else
      facet_ptr->outwards_unit_normal = candidate_normal;
  }

  // Create frontiers given group black vertices
  for (int i = 0; i < num_groups; i++) {
    vector<FacetPtr> group_facets =
        findGroupFacetsFromVertices(node, bv_groups.at(i));
    if (group_facets.empty())
      continue;

    findNbhdFacets(group_facets);
    vector<vector<FacetPtr>> linked_groups;
    for (FacetPtr facet : group_facets) {
      if (facet->linked)
        continue;

      vector<FacetPtr> linked_group_facets;
      deque<FacetPtr> pending_facets;
      pending_facets.push_back(facet);
      while (!pending_facets.empty()) {
        FacetPtr current_facet = pending_facets.front();
        pending_facets.pop_front();
        if (current_facet->linked)
          continue;
        linked_group_facets.push_back(current_facet);
        current_facet->linked = true;
        for (FacetPtr nbhd : current_facet->nbhd_facets) {
          if (!nbhd->linked)
            pending_facets.push_back(nbhd);
        }
      }

      linked_groups.push_back(linked_group_facets);
    }

    int num_linked_groups = linked_groups.size();
    for (int j = 0; j < num_linked_groups; j++) {
      vector<FrontierPtr> frontiers = splitFrontier(node, linked_groups.at(j));
      // vector<FrontierPtr> frontiers = splitFrontier(node, group_facets);
      for (FrontierPtr f : frontiers) {
        // ROS_ERROR("Pushed a new frontier, proj_center(%f, %f, %f)",
        // f->proj_center(0),
        //           f->proj_center(1), f->proj_center(2));
        node->frontiers.push_back(f);
        // for (auto f : pending_frontiers) {
        //   ROS_INFO_STREAM("Frontier pos: " << f->next_node_pos.transpose());
        // }
      }
    }
  }

  // Add bv_group for those connecting bvs having big diff in dis_to_center
  vector<FacetPtr> jump_facets;
  for (FacetPtr facet : node->facets) {
    if (facet->valid)
      continue;
    if (facetOnCeilOrFloor(facet))
      continue;

    int count_jump = 0;
    // bool bad = false;
    for (int i = 0; i < 3; i++) {
      VertexPtr v1 = facet->vertices.at(i);
      // if(onCeilOrFloor(v1->coord) != 0){
      //   bad = true;
      //   break;
      // }
      VertexPtr v2 = facet->vertices.at((i + 1) % 3);
      if (v1->dis_to_center > _frontier_jump_threshold * v2->dis_to_center ||
          v2->dis_to_center > _frontier_jump_threshold * v1->dis_to_center) {
        count_jump++;
      }
    }
    // if (bad) continue;
    if (count_jump > 1) {
      jump_facets.push_back(facet);
      facet->valid = true;
    }
  }

  findNbhdFacets(jump_facets);
  vector<vector<FacetPtr>> linked_groups;
  for (FacetPtr facet : jump_facets) {
    if (facet->linked)
      continue;

    vector<FacetPtr> linked_group_facets;
    deque<FacetPtr> pending_facets;
    pending_facets.push_back(facet);
    while (!pending_facets.empty()) {
      FacetPtr current_facet = pending_facets.front();
      pending_facets.pop_front();
      if (current_facet->linked)
        continue;
      linked_group_facets.push_back(current_facet);
      current_facet->linked = true;
      for (FacetPtr nbhd : current_facet->nbhd_facets) {
        if (!nbhd->linked)
          pending_facets.push_back(nbhd);
      }
    }

    linked_groups.push_back(linked_group_facets);
  }

  int num_linked_groups = linked_groups.size();
  for (int j = 0; j < num_linked_groups; j++) {
    FrontierPtr new_fron = new Frontier(linked_groups.at(j), node);
    if (initFrontier(new_fron))
      node->frontiers.push_back(new_fron);
    // ROS_INFO("Added a jump frontier!");
    // for (auto f : pending_frontiers) {
    //   ROS_INFO_STREAM("Frontier pos: " << f->next_node_pos.transpose());
    // }
  }

  // Wish to expand frontier with more facets first
  sort(node->frontiers.begin(), node->frontiers.end(), compareFrontier);
  for (FrontierPtr f : node->frontiers) {
    // ROS_INFO("This frontier size: %d", f->facets.size());
    if (f->facets.empty())
      continue;
    pending_frontiers.push_back(f);
    int ind = loop_candidate_frontiers.size();
    f->index = ind;
    loop_candidate_frontiers.push_back(f);
  }
}

vector<FrontierPtr>
SkeletonFinder::splitFrontier(NodePtr node, vector<FacetPtr> group_facets) {
  vector<FrontierPtr> frontiers;
  // bool test_node = false;
  // Eigen::Vector3d position(-20.912447, 5.552020, 1.432291);
  // if (getDis(position, node->coord) < 0.5) {
  //   ROS_ERROR("Initing testing node!");
  //   getchar();
  //   test_node = true;
  // }
  if ((int)group_facets.size() <= _max_facets_grouped) {
    // if (!test_node) {
    // FrontierPtr new_frontier = new Frontier(group_facets, node);
    // frontiers.push_back(new_frontier);
    // if (!initFrontier(new_frontier)) {
    //   ROS_ERROR("1. Init frontier failed!!! Change to use alternative
    //   proj_center");
    // }
    // }
    // else {
    Eigen::Vector3d avg_normal = Eigen::Vector3d::Zero();
    vector<FacetPtr> filtered_group_facets;
    for (FacetPtr f : group_facets) {
      avg_normal += f->outwards_unit_normal;
    }
    avg_normal.normalize();
    for (FacetPtr f : group_facets) {
      if (f->nbhd_facets.size() < 2) {
        double angle = acos(avg_normal.dot(f->outwards_unit_normal));
        if (angle > M_PI / 2.5) {
          vector<FacetPtr> single_facet;
          single_facet.push_back(f);
          FrontierPtr new_frontier = new Frontier(single_facet, node);
          frontiers.push_back(new_frontier);
          if (!initFrontier(new_frontier)) {
            // ROS_ERROR("1-2. Init frontier failed!!! Change to use alternative
            // center and normal");
          }
          continue;
        }
      }
      filtered_group_facets.push_back(f);
    }
    FrontierPtr new_frontier = new Frontier(filtered_group_facets, node);
    frontiers.push_back(new_frontier);
    if (!initFrontier(new_frontier)) {
      // ROS_ERROR("1-2. Init frontier failed!!! Change to use alternative
      // proj_center");
    }
    // }
  } else {
    // findNbhdFacets(group_facets);
    for (FacetPtr facet : group_facets) {
      if (facet->visited)
        continue;

      facet->visited = true;
      Eigen::Vector3d normal = Eigen::Vector3d::Zero();
      vector<FacetPtr> small_group_facets;
      deque<FacetPtr> pending_facets;
      pending_facets.push_back(facet);

      while (!pending_facets.empty() &&
             (int)small_group_facets.size() < _max_facets_grouped) {
        // while (!pending_facets.empty()) {
        FacetPtr f = pending_facets.front();
        pending_facets.pop_front();

        if (!small_group_facets.empty()) {
          if (acos(f->outwards_unit_normal.dot(normal)) <
              M_PI / _frontier_split_threshold) {
            normal =
                (normal * small_group_facets.size() + f->outwards_unit_normal) /
                (small_group_facets.size() + 1);
            normal.normalize();
          } else
            continue;
        } else {
          normal = f->outwards_unit_normal;
        }

        /*
        bool fit = true;
        if (!small_group_facets.empty()) {
          for (FacetPtr f_in_group : small_group_facets) {
            if
        (acos(f->outwards_unit_normal.dot(f_in_group->outwards_unit_normal)) >
        M_PI / 3.0
        * 2.0) { fit = false; break;
            }
          }
        }
        if (!fit) continue;
        */

        f->visited = true;
        small_group_facets.push_back(f);
        for (FacetPtr f_nbhd : f->nbhd_facets) {
          if (f_nbhd->visited)
            continue;
          pending_facets.push_back(f_nbhd);
        }
      }

      FrontierPtr new_frontier = new Frontier(small_group_facets, node);
      frontiers.push_back(new_frontier);
      if (!initFrontier(new_frontier)) {
        // ROS_ERROR("2. Init frontier failed!!! Change to use alternative
        // center and normal");
      }
    }
  }
  // if (frontiers.size() > 1) ROS_WARN("Split frontier into %d frontiers",
  // frontiers.size());
  return frontiers;
}

void SkeletonFinder::findNbhdFacets(vector<FacetPtr> facets) {
  int num_facet = facets.size();
  for (int i = 0; i < num_facet - 1; i++) {
    for (int j = i + 1; j < num_facet; j++) {
      FacetPtr f1 = facets.at(i);
      FacetPtr f2 = facets.at(j);

      if (f1->nbhd_facets.size() == 3)
        break;

      int same_vertices_count = 0;
      for (int s = 0; s < 3; s++) {
        for (int t = 0; t < 3; t++) {
          if (isSamePos(f1->vertices.at(s)->coord, f2->vertices.at(t)->coord)) {
            same_vertices_count++;
            break;
          }
        }
      }
      if (same_vertices_count == 2) {
        // if (same_vertices_count >= 1) {
        f1->nbhd_facets.push_back(f2);
        f2->nbhd_facets.push_back(f1);
      }
    }
  }
}

vector<FacetPtr>
SkeletonFinder::findGroupFacetsFromVertices(NodePtr node,
                                            vector<VertexPtr> group_bv) {
  vector<FacetPtr> group_facets;
  for (FacetPtr f : node->facets) {
    // int num_bv_included = 0;
    bool good = true;
    for (VertexPtr v_facet : f->vertices) {
      bool notIncluded = true;
      for (VertexPtr v_group : group_bv) {
        if (!v_group->critical)
          continue;
        if (isSamePos(v_facet->coord, v_group->coord)) {
          notIncluded = false;
          break;
        }
      }
      if (notIncluded) {
        good = false;
        break;
      }
      // if (!notIncluded) {
      //   num_bv_included++;
      // }
    }
    // All of the three vertices are included in the group_bv
    if (good) {
      // if (num_bv_included >= 3) {
      // Exclude facets completely on ceil or floor to improve frontier accuracy
      if ((onCeilOrFloor(f->vertices.at(0)->coord) == 1 &&
           onCeilOrFloor(f->vertices.at(1)->coord) == 1 &&
           onCeilOrFloor(f->vertices.at(2)->coord) == 1) ||
          (onCeilOrFloor(f->vertices.at(0)->coord) == -1 &&
           onCeilOrFloor(f->vertices.at(1)->coord) == -1 &&
           onCeilOrFloor(f->vertices.at(2)->coord) == -1))
        continue;
      group_facets.push_back(f);
      f->valid = true;
    }
  }
  return group_facets;
}

void SkeletonFinder::findFlowBack(NodePtr node) {
  if (node->seed_frontier == NULL)
    return;
  // vector<int> loop_frontier_counter(loop_candidate_frontiers.size(), 0);
  int size = loop_candidate_frontiers.size();
  vector<vector<Eigen::Vector3d>> flow_back_frontier_log(size);
  vector<FrontierPtr> frontiers(size);
  vector<int> pending_frontier_index;
  vector<int> collision_node_index_log;

  // Count number of contact black vertices on each frontiers
  for (VertexPtr v : node->black_vertices) {
    // only process vertices collide with other polyhedrons
    if (v->collision_node_index < 0)
      continue;
    // hit_on_pcl is on seed_frontier which is already connected
    if (checkPtOnFrontier(node->seed_frontier, v->coord)) {
      // ROS_WARN("hit_on_pcl is on seed_frontier");
      continue;
    }

    FrontierPtr loop_ftr =
        findFlowBackFrontier(v->coord, v->collision_node_index);
    if (loop_ftr == NULL) {
      // ROS_ERROR("Cannot find a loop frontier!!!");
      collision_node_index_log.push_back(v->collision_node_index);
      continue;
    }
    // loop_frontier_counter.at(loop_ftr->index)++;
    flow_back_frontier_log.at(loop_ftr->index).push_back(v->coord);
    frontiers.at(loop_ftr->index) = loop_ftr;
  }

  // Sort decreasingly: flowback frontiers with more hits first
  for (int i = 0; i < size; i++) {
    if (flow_back_frontier_log.at(i).empty())
      continue;
    if ((int)flow_back_frontier_log.at(i).size() <
        _min_flowback_creation_threshold) {
      // ROS_ERROR("Flowback size: %d", flow_back_frontier_log.at(i).size());
      // ROS_ERROR("Radius: %f",
      // getVerticesRadius(flow_back_frontier_log.at(i)));
      if (getVerticesRadius(flow_back_frontier_log.at(i)) <
          _min_flowback_creation_radius_threshold)
        continue;
    }
    if (pending_frontier_index.empty())
      pending_frontier_index.push_back(i);
    else {
      vector<int>::iterator it;
      for (it = pending_frontier_index.begin();
           it != pending_frontier_index.end(); it++) {
        if (flow_back_frontier_log.at(*it).size() <
            flow_back_frontier_log.at(i).size()) {
          pending_frontier_index.insert(it, i);
          break;
        }
        pending_frontier_index.push_back(i);
        break;
      }
    }
  }

  // Start flowback
  int size_pending_frontier = pending_frontier_index.size();
  // ROS_ERROR("size_pending_frontier: %d", size_pending_frontier);
  vector<Eigen::Vector3d> connected_node_pos;
  if (node->seed_frontier->gate_node->connected_Node_ptr.empty()) {
    connected_node_pos.push_back(node->seed_frontier->master_node->coord);
    // ROS_INFO("gate node con empty: push_back a pos");
  } else {
    for (NodePtr gate_con_node :
         node->seed_frontier->gate_node->connected_Node_ptr) {
      connected_node_pos.push_back(gate_con_node->coord);
      // ROS_INFO("gate node con nonempty: push_back a pos");
    }
  }
  for (int i = 0; i < size_pending_frontier; i++) {
    // ROS_ERROR("%d ,", pending_frontier_index.at(i));
    int ind = pending_frontier_index.at(i);
    FrontierPtr flowback_frontier = frontiers.at(ind);

    // Unsafe loop: loop seg is not obstacle-free
    Eigen::Vector3d end_pt_on_frontier;
    if (flowback_frontier->gate_node == NULL)
      end_pt_on_frontier = flowback_frontier->proj_center;
    else
      end_pt_on_frontier = flowback_frontier->gate_node->coord;
    if (!checkSegClear(node->coord, end_pt_on_frontier) ||
        !checkSegClear(flowback_frontier->master_node->coord,
                       end_pt_on_frontier)) {
      // ROS_ERROR("Flowback seg not clear!");
      continue;
    }

    // Bad loop: loop only contains 4 nodes
    if (_bad_loop) {
      bool bad_loop = false;
      for (Eigen::Vector3d pos : connected_node_pos) {
        // ROS_INFO("Checking connected node pos: (%f, %f, %f)", pos(0), pos(1),
        // pos(2));
        if (flowback_frontier->gate_node == NULL ||
            flowback_frontier->gate_node->rollbacked) {
          // ROS_INFO("flowback_frontier doesn't have gate node");
          if (isSamePos(pos, flowback_frontier->master_node->coord)) {
            bad_loop = true;
            // ROS_INFO("isSamePos");
          }
          // ROS_INFO("NOT isSamePos");
        } else {
          // ROS_INFO("flowback_frontier has gate node");
          for (NodePtr frontier_con_node :
               flowback_frontier->gate_node->connected_Node_ptr) {
            if (isSamePos(pos, frontier_con_node->coord)) {
              bad_loop = true;
              // ROS_INFO("isSamePos");
              break;
            }
            // ROS_INFO("NOT isSamePos");
          }
        }
        if (bad_loop)
          break;
      }
      if (bad_loop) {
        // ROS_ERROR("BAD LOOP!!!");
        continue;
      }
    }

    // Create flowback: new gate node if necessary
    if (flowback_frontier->gate_node == NULL) {
      NodePtr new_gate =
          new Node(flowback_frontier->proj_center, flowback_frontier, true);
      if (initNode(new_gate)) {
        flowback_frontier->gate_node = new_gate;
        flowback_frontier->master_node->connected_Node_ptr.push_back(new_gate);
        new_gate->connected_Node_ptr.push_back(flowback_frontier->master_node);
        node->connected_Node_ptr.push_back(new_gate);
        new_gate->connected_Node_ptr.push_back(node);
      } else
        continue;
      // flowback_frontier->gate_node = new_gate;
      // flowback_frontier->master_node->connected_Node_ptr.push_back(new_gate);
      // new_gate->connected_Node_ptr.push_back(flowback_frontier->master_node);
      // initNode(new_gate);
      // node->connected_Node_ptr.push_back(new_gate);
      // new_gate->connected_Node_ptr.push_back(node);
    } else if (flowback_frontier->gate_node->rollbacked) {
      flowback_frontier->gate_node->rollbacked = false;
      recordNode(flowback_frontier->gate_node);
      flowback_frontier->master_node->connected_Node_ptr.push_back(
          flowback_frontier->gate_node);
      flowback_frontier->gate_node->connected_Node_ptr.push_back(
          flowback_frontier->master_node);
      node->connected_Node_ptr.push_back(flowback_frontier->gate_node);
      flowback_frontier->gate_node->connected_Node_ptr.push_back(node);
    } else {
      node->connected_Node_ptr.push_back(flowback_frontier->gate_node);
      flowback_frontier->gate_node->connected_Node_ptr.push_back(node);
    }

    for (NodePtr frontier_con_node :
         flowback_frontier->gate_node->connected_Node_ptr) {
      connected_node_pos.push_back(frontier_con_node->coord);
    }
  }
}

bool SkeletonFinder::initFrontier(FrontierPtr frontier) {
  // Set proj_center
  bool proj_center_found = false;

  // Line equation:
  // x = x0 + t * nx
  // y = y0 + t * ny
  // z = z0 + t * nz
  double x0 = frontier->avg_center(0);
  double y0 = frontier->avg_center(1);
  double z0 = frontier->avg_center(2);
  double nx = frontier->outwards_unit_normal(0);
  double ny = frontier->outwards_unit_normal(1);
  double nz = frontier->outwards_unit_normal(2);

  int num_facet = frontier->facets.size();
  for (int i = 0; i < num_facet; i++) {
    double a = frontier->facets.at(i)->plane_equation(0);
    double b = frontier->facets.at(i)->plane_equation(1);
    double c = frontier->facets.at(i)->plane_equation(2);
    double d = frontier->facets.at(i)->plane_equation(3);

    double t = -(a * x0 + b * y0 + c * z0 + d) / (a * nx + b * ny + c * nz);

    Eigen::Vector3d intersection =
        frontier->avg_center + t * frontier->outwards_unit_normal;

    Eigen::Vector3d coord1 = frontier->facets.at(i)->vertices.at(0)->coord;
    Eigen::Vector3d coord2 = frontier->facets.at(i)->vertices.at(1)->coord;
    Eigen::Vector3d coord3 = frontier->facets.at(i)->vertices.at(2)->coord;

    Eigen::Vector3d cross1 = (coord2 - coord1).cross(intersection - coord1);
    Eigen::Vector3d cross2 = (coord3 - coord2).cross(intersection - coord2);
    Eigen::Vector3d cross3 = (coord1 - coord3).cross(intersection - coord3);
    if (cross1(0) * cross2(0) > 0 && cross2(0) * cross3(0) > 0 &&
        cross3(0) * cross1(0) > 0) {
      frontier->proj_center = intersection;
      // frontier->proj_facet_normal =
      // frontier->facets.at(i)->outwards_unit_normal;
      // frontier->proj_facet_center = frontier->facets.at(i)->center;
      frontier->proj_facet = frontier->facets.at(i);
      Eigen::Vector3d normal1 = frontier->outwards_unit_normal;
      Eigen::Vector3d normal2 = frontier->facets.at(i)->outwards_unit_normal;
      frontier->cos_theta =
          normal1.dot(normal2) / (normal1.norm() * normal2.norm());
      proj_center_found = true;
      break;
    }
  }

  if (!proj_center_found) {
    // ROS_ERROR("size of frontier facets: %d", frontier->facets.size());
    double min_angle = M_PI;
    FacetPtr best_facet;
    for (FacetPtr f : frontier->facets) {
      double angle =
          acos(frontier->outwards_unit_normal.dot(f->outwards_unit_normal));
      if (angle < min_angle) {
        min_angle = angle;
        best_facet = f;
      }
    }
    frontier->proj_facet = best_facet;
    frontier->proj_center = best_facet->center;
    // frontier->proj_facet_normal = best_facet->outwards_unit_normal;
    // frontier->proj_facet_center = best_facet->center;
  }

  // Set vertices
  for (FacetPtr facet : frontier->facets) {
    for (VertexPtr v_facet : facet->vertices) {
      bool exist = false;
      if (!frontier->vertices.empty()) {
        for (VertexPtr v_frontier : frontier->vertices) {
          if (isSamePos(v_facet->coord, v_frontier->coord)) {
            exist = true;
            break;
          }
        }
      }
      if (!exist)
        frontier->vertices.push_back(v_facet);
    }
  }

  return proj_center_found;
}

bool SkeletonFinder::checkPtInPolyhedron(NodePtr node, Eigen::Vector3d pt) {
  int intersection_count = 0;
  double x = pt(0);
  double y = pt(1);
  double z = pt(2);

  int num_ftr = node->facets.size();
  for (int i = 0; i < num_ftr; i++) {
    FacetPtr inter_ftr_ptr = node->facets.at(i);

    Eigen::Vector3d coord1 = inter_ftr_ptr->vertices.at(0)->coord;
    Eigen::Vector3d coord2 = inter_ftr_ptr->vertices.at(1)->coord;
    Eigen::Vector3d coord3 = inter_ftr_ptr->vertices.at(2)->coord;

    double a = inter_ftr_ptr->plane_equation(0);
    double b = inter_ftr_ptr->plane_equation(1);
    double c = inter_ftr_ptr->plane_equation(2);
    double d = inter_ftr_ptr->plane_equation(3);

    double inter_z = -(a * x + b * y + d) / c;
    Eigen::Vector3d pt_to_judge;
    pt_to_judge << x, y, inter_z;

    if (inter_z >= z) {
      Eigen::Vector3d cross1 = (coord2 - coord1).cross(pt_to_judge - coord1);
      Eigen::Vector3d cross2 = (coord3 - coord2).cross(pt_to_judge - coord2);
      Eigen::Vector3d cross3 = (coord1 - coord3).cross(pt_to_judge - coord3);
      if (cross1(0) * cross2(0) > 0 && cross2(0) * cross3(0) > 0 &&
          cross3(0) * cross1(0) > 0) {
        intersection_count++;
      }
    }
  }

  if (intersection_count % 2 == 0)
    return false;
  else
    return true;
}

void SkeletonFinder::verifyFrontier(FrontierPtr ftr_ptr) {
  // ROS_ERROR("Start verifyFrontier");

  Eigen::Vector3d raycast_start_pt =
      ftr_ptr->proj_center + 2.0 * ftr_ptr->outwards_unit_normal *
                                 _search_margin; // / ftr_ptr->cos_theta;
  // Eigen::Vector3d raycast_start_pt =
  //     ftr_ptr->proj_center + ftr_ptr->outwards_unit_normal * _search_margin;
  // ROS_ERROR("raycast_start_pt(%f, %f, %f)", raycast_start_pt(0),
  // raycast_start_pt(1),
  //           raycast_start_pt(2));

  pair<double, int> rs_result = radiusSearch(raycast_start_pt);
  if (rs_result.first < _search_margin) {
    ftr_ptr->valid = false;
    return;
  }
  // ROS_ERROR("raycast_start_pt OK dis to obstable: %f",
  // radiusSearch(raycast_start_pt).first);

  pair<Vector3d, int> raycast_result =
      raycast(raycast_start_pt, ftr_ptr->outwards_unit_normal,
              _max_expansion_ray_length);
  Eigen::Vector3d hit_on_pcl = raycast_result.first;

  // raycast into a long corridor
  if (hit_on_pcl == raycast_start_pt) {
    // ROS_INFO("verifyFrontiers: raycast into a long corridor");
    Eigen::Vector3d new_node_candidate =
        ftr_ptr->proj_center +
        0.5 * _max_expansion_ray_length * ftr_ptr->outwards_unit_normal;
    if (checkWithinBbx(new_node_candidate)) {
      ftr_ptr->next_node_pos = new_node_candidate;
      ftr_ptr->valid = true;
      // ROS_INFO("verifyFrontiers: and set new node at max limit");
    } else {
      // ROS_INFO("verifyFrontiers: new node candidate exceeds bbx");
    }
  }
  // normal case
  else if (getDis(hit_on_pcl, ftr_ptr->proj_center) >
           _frontier_creation_threshold) {
    ftr_ptr->valid = true;
    ftr_ptr->next_node_pos = (hit_on_pcl + ftr_ptr->proj_center) / 2;
    // ROS_INFO("verifyFrontiers: valid. New node (%f, %f, %f)",
    // ftr_ptr->next_node_pos(0),
    //          ftr_ptr->next_node_pos(1), ftr_ptr->next_node_pos(2));
  } else {
    // ROS_INFO("verifyFrontiers: hit_on_pcl too near to frontier");
  }
}

void SkeletonFinder::addFacetsToPcl(NodePtr nodePtr) {
  pcl::PointCloud<pcl::PointXYZ> poly_pcl;
  int num_facet = nodePtr->facets.size();
  for (int i = 0; i < num_facet; i++) {
    FacetPtr facet = nodePtr->facets.at(i);

    vector<Eigen::Vector3d> start_list;
    vector<double> length_list;

    Eigen::Vector3d top_vertex = facet->vertices.at(0)->coord;
    Eigen::Vector3d left_vertex = facet->vertices.at(1)->coord;
    Eigen::Vector3d right_vertex = facet->vertices.at(2)->coord;

    Eigen::Vector3d left_to_top = top_vertex - left_vertex;
    Eigen::Vector3d left_to_top_unit = left_to_top / left_to_top.norm();
    Eigen::Vector3d left_to_right = right_vertex - left_vertex;
    Eigen::Vector3d right_to_top = top_vertex - right_vertex;
    Eigen::Vector3d right_to_top_unit = right_to_top / right_to_top.norm();
    Eigen::Vector3d right_to_left = left_vertex - right_vertex;

    double theta = acos(left_to_top.dot(left_to_right) /
                        (left_to_top.norm() * left_to_right.norm()));
    double theta_right = acos(right_to_top.dot(right_to_left) /
                              (right_to_top.norm() * right_to_left.norm()));
    double step_length = _resolution / 2;
    // double step_length = _search_margin / 2;
    double start_pt_step = step_length / sin(theta);
    double end_pt_step = step_length / sin(theta_right);

    int num_start_pt = ceil(left_to_top.norm() / start_pt_step);
    double length;
    for (int j = 0; j < num_start_pt; j++) {
      Eigen::Vector3d start =
          left_vertex + j * start_pt_step * left_to_top_unit;
      Eigen::Vector3d end = right_vertex + j * end_pt_step * right_to_top_unit;
      start_list.push_back(start);
      length = (end - start).norm();
      length_list.push_back(length);
    }

    Eigen::Vector3d direction_unit = right_vertex - left_vertex;
    direction_unit.normalize();
    for (int j = 0; j < num_start_pt; j++) {
      int num_pts = ceil(length_list.at(j) / step_length);
      for (int k = 0; k < num_pts; k++) {
        Eigen::Vector3d pt_to_push =
            start_list.at(j) + k * direction_unit * step_length;
        for (int s = 0; s < 4; s++) {
          // map_pcl.points.push_back(pcl::PointXYZ(pt_to_push(0),
          // pt_to_push(1), pt_to_push(2)));
          poly_pcl.points.push_back(
              pcl::PointXYZ(pt_to_push(0), pt_to_push(1), pt_to_push(2)));
          pt_to_push += (-facet->outwards_unit_normal) * step_length;
        }
      }
      Eigen::Vector3d end_to_push =
          start_list.at(j) + length_list.at(j) * direction_unit;
      for (int s = 0; s < 4; s++) {
        // map_pcl.points.push_back(pcl::PointXYZ(end_to_push(0),
        // end_to_push(1), end_to_push(2)));
        poly_pcl.points.push_back(
            pcl::PointXYZ(end_to_push(0), end_to_push(1), end_to_push(2)));
        end_to_push += (-facet->outwards_unit_normal) * step_length;
      }
    }
    // map_pcl.points.push_back(pcl::PointXYZ(top_vertex(0), top_vertex(1),
    // top_vertex(2)));
    poly_pcl.points.push_back(
        pcl::PointXYZ(top_vertex(0), top_vertex(1), top_vertex(2)));
  }

  shared_ptr<pcl::search::KdTree<pcl::PointXYZ>> new_poly(
      new pcl::search::KdTree<pcl::PointXYZ>);
  // new_poly.reset(new pcl::search::KdTree<pcl::PointXYZ>);
  new_poly->setInputCloud(poly_pcl.makeShared());
  kdtreesForPolys.push_back(new_poly);

  // ROS_WARN_STREAM("map pcl size: " << map_pcl.points.size());

  /*
  kdtreeForMapPtr.reset(new pcl::search::KdTree<pcl::PointXYZ>);
  ros::Time begin = ros::Time::now();
  // begin = ros::Time::now();
  kdtreeForMapPtr->setInputCloud(map_pcl.makeShared());
  ros::Time finish = ros::Time::now();
  // finish = ros::Time::now();
  ROS_ERROR("Updating map kdtree takes time: %f", (finish - begin).toSec());
  kdtree_timing += (finish - begin).toSec();
  */
}

bool SkeletonFinder::processFrontier(FrontierPtr curFtrPtr) {
  // ROS_ERROR("Start process frontier");
  // double direction_theta =
  //     atan2(curFtrPtr->outwards_unit_normal(1),
  //     curFtrPtr->outwards_unit_normal(0)) * 180 / M_PI;

  // Gate node: midpoint of frontier
  NodePtr gate;
  if (curFtrPtr->gate_node == NULL) {
    gate = new Node(curFtrPtr->proj_center, curFtrPtr, true);
    curFtrPtr->gate_node = gate;
  } else {
    gate = curFtrPtr->gate_node;
  }

  bool floor = checkFloor(gate);
  bool bbx = checkWithinBbx(gate->coord);
  if (!floor || !bbx) {
    gate->rollbacked = true;
    // if (!floor) {
    //   ROS_INFO("processFrontier: no floor");
    // }
    // if (!bbx) {
    //   ROS_INFO("processFrontier: outside bbx");
    // }
    // ROS_INFO("processFrontier: gate coord: (%f, %f, %f)", gate->coord(0),
    // gate->coord(1),
    //          gate->coord(2));
    return false;
  }

  // Center node
  NodePtr new_node = new Node(curFtrPtr->next_node_pos, curFtrPtr);
  bool init_success = initNode(new_node);

  if (init_success) {
    initNode(gate);
    // double-sided pointers
    curFtrPtr->master_node->connected_Node_ptr.push_back(gate);
    gate->connected_Node_ptr.push_back(curFtrPtr->master_node);
    gate->connected_Node_ptr.push_back(new_node);
    new_node->connected_Node_ptr.push_back(gate);
  } else {
    // ROS_INFO("processFrontier: new node init fails");
    new_node->rollbacked = true;
    if (gate->connected_Node_ptr.empty()) {
      gate->rollbacked = true;
      curFtrPtr->valid = false;
      for (auto f : curFtrPtr->facets) {
        f->valid = false;
      }
    }
  }
  // ROS_WARN("Number of connected node: %d",
  // new_node->connected_Node_ptr.size()); for (NodePtr con :
  // new_node->connected_Node_ptr) {
  //   ROS_INFO("Connected node pos (%f, %f, %f)", con->coord(0), con->coord(1),
  //   con->coord(2));
  // }

  return init_success;
}

vector<Eigen::Vector3d> SkeletonFinder::findPath(Eigen::Vector3d start,
                                                 Eigen::Vector3d target) {
  Eigen::Vector3d start_astar = Eigen::Vector3d::Zero();
  Eigen::Vector3d target_astar = Eigen::Vector3d::Zero();
  vector<Eigen::Vector3d> path;

  // for (NodePtr node : NodeList) {
  //   if (node->rollbacked) continue;
  //   nodes_pcl.points.push_back(pcl::PointXYZ(node->coord(0), node->coord(1),
  //   node->coord(2)));
  // }
  kdtreeForNodes.setInputCloud(nodes_pcl.makeShared());

  pcl::PointXYZ pcl_start(start(0), start(1), start(2));
  pointIdxRadiusSearchForNodes.clear();
  pointRadiusSquaredDistanceForNodes.clear();
  kdtreeForNodes.nearestKSearch(pcl_start, 5, pointIdxRadiusSearchForNodes,
                                pointRadiusSquaredDistanceForNodes);
  for (std::size_t i = 0; i < pointIdxRadiusSearchForNodes.size(); ++i) {
    Eigen::Vector3d nbhd;
    nbhd << nodes_pcl[pointIdxRadiusSearchForNodes[i]].x,
        nodes_pcl[pointIdxRadiusSearchForNodes[i]].y,
        nodes_pcl[pointIdxRadiusSearchForNodes[i]].z;
    if (checkSegClear(start, nbhd)) {
      start_astar = nbhd;
      break;
    }
  }

  pcl::PointXYZ pcl_target(target(0), target(1), target(2));
  pointIdxRadiusSearchForNodes.clear();
  pointRadiusSquaredDistanceForNodes.clear();
  kdtreeForNodes.nearestKSearch(pcl_target, 10, pointIdxRadiusSearchForNodes,
                                pointRadiusSquaredDistanceForNodes);
  for (std::size_t i = 0; i < pointIdxRadiusSearchForNodes.size(); ++i) {
    Eigen::Vector3d nbhd;
    nbhd << nodes_pcl[pointIdxRadiusSearchForNodes[i]].x,
        nodes_pcl[pointIdxRadiusSearchForNodes[i]].y,
        nodes_pcl[pointIdxRadiusSearchForNodes[i]].z;
    if (checkSegClear(target, nbhd)) {
      target_astar = nbhd;
      break;
    }
  }

  if (start_astar == Eigen::Vector3d::Zero() ||
      target_astar == Eigen::Vector3d::Zero()) {
    ROS_ERROR("Can't find nodes on skeleton to connect!");
    return path;
  }

  vector<Eigen::Vector3d> astar_path =
      findPathByAStar(start_astar, target_astar);
  if (!astar_path.empty()) {
    path.push_back(start);
    for (Eigen::Vector3d waypoint : astar_path)
      path.push_back(waypoint);
    path.push_back(target);
    ROS_INFO("Path found.");
    path_finder.visNodes();
    path_finder.visConnections();
    return path;
  } else {
    ROS_ERROR("Path not found!");
    return path;
  }
}

vector<Eigen::Vector3d>
SkeletonFinder::findPathByAStar(Eigen::Vector3d start, Eigen::Vector3d target) {
  vector<a_star::NodePtr> as_nodes;
  for (NodePtr node : NodeList) {
    if (node->rollbacked)
      continue;

    vector<Eigen::Vector3d> connected_node_pos;
    for (NodePtr con_node : node->connected_Node_ptr) {
      // if (con_node->rollbacked) continue;
      connected_node_pos.push_back(con_node->coord);
    }
    a_star::NodePtr as_node = new a_star::Node(node->coord, connected_node_pos);
    as_nodes.push_back(as_node);
  }

  path_finder.init(as_nodes);
  int result = path_finder.search(start, target);
  vector<Eigen::Vector3d> path;
  if (result == 1) { // REACH_END
    // ROS_INFO("A_star: REACH_END!");
    path = path_finder.getPath();
  } else {
    // ROS_ERROR("A_star: NO_PATH!");
  }
  return path;
}

/* ------------------ Utility functions ----------------- */
bool SkeletonFinder::isSamePos(Eigen::Vector3d pos1, Eigen::Vector3d pos2) {
  return getDis(pos1, pos2) < 1e-4;
}

VertexPtr SkeletonFinder::getVertexFromDire(NodePtr node,
                                            Eigen::Vector3d dire) {
  for (VertexPtr v : node->black_vertices) {
    if (isSamePos(v->dire_unit_sphere, dire))
      return v;
  }
  for (VertexPtr v : node->white_vertices) {
    if (isSamePos(v->dire_unit_sphere, dire))
      return v;
  }
  return NULL;
}

bool SkeletonFinder::checkFloor(NodePtr node) {
  Eigen::Vector3d downwards(0, 0, -1);
  pair<Vector3d, int> raycast_result =
      raycastOnRawMap(node->coord, downwards, _max_ray_length);
  // ROS_ERROR("raycast result: %d", raycast_result.second);
  // ROS_ERROR("raycast result: (%f, %f, %f)", raycast_result.first(0),
  // raycast_result.first(1),
  //           raycast_result.first(2));
  if (raycast_result.second == -2) {
    // ROS_INFO("Floor not found within max_ray_length");
    return false;
  }
  double floor_height = raycast_result.first(2);

  // First node case
  if (node->seed_frontier == NULL) {
    node->dis_to_floor = floor_height;
    return true;
  }
  if (!node->isGate) {
    Eigen::Vector3d mid = (node->coord + node->seed_frontier->proj_center) / 2;
    pair<Vector3d, int> raycast_result_mid =
        raycastOnRawMap(mid, downwards, _max_ray_length);
    if (raycast_result_mid.second == -2) {
      // ROS_INFO("Floor not found within max_ray_length");
      return false;
    }

    Eigen::Vector3d mid2 = (node->seed_frontier->proj_center +
                            node->seed_frontier->master_node->coord) /
                           2;
    pair<Vector3d, int> raycast_result_mid2 =
        raycastOnRawMap(mid2, downwards, _max_ray_length);
    if (raycast_result_mid2.second == -2) {
      // ROS_INFO("Floor not found within max_ray_length");
      return false;
    }
  }

  double parent_floor_height = node->seed_frontier->master_node->dis_to_floor;
  if (fabs(floor_height - parent_floor_height) > _max_height_diff) {
    // ROS_INFO("large diff from parent");
    // ROS_INFO("parent_floor_height: %f", parent_floor_height);
    // ROS_INFO("this floor_height: %f", floor_height);
    return false;
  }

  node->dis_to_floor = floor_height;

  return true;
}

int SkeletonFinder::onCeilOrFloor(Eigen::Vector3d p) {
  // On ceil
  if (fabs(p(2) - _z_max) < _search_margin)
    return 1;
  // On floor
  if (fabs(p(2) - _z_min) < _search_margin)
    return -1;
  // Not on ceil or floor
  return 0;
}

bool SkeletonFinder::facetOnCeilOrFloor(FacetPtr f) {
  return (onCeilOrFloor(f->vertices.at(0)->coord) == -1 &&
          onCeilOrFloor(f->vertices.at(1)->coord) == -1 &&
          onCeilOrFloor(f->vertices.at(2)->coord) == -1) ||
         (onCeilOrFloor(f->vertices.at(0)->coord) == 1 &&
          onCeilOrFloor(f->vertices.at(1)->coord) == 1 &&
          onCeilOrFloor(f->vertices.at(2)->coord) == 1);
}

bool SkeletonFinder::checkSegClear(Eigen::Vector3d pos1, Eigen::Vector3d pos2) {
  double length = (pos2 - pos1).norm();
  double step_length = _resolution;

  Eigen::Vector3d step = step_length * (pos2 - pos1) / length;
  int num_steps = ceil(length / step_length);

  Eigen::Vector3d begin_pos = pos1;
  for (int i = 0; i < num_steps; i++) {
    Eigen::Vector3d check_pos = begin_pos + i * step;
    if (collisionCheck(check_pos, _search_margin)) {
      return false;
    }
  }
  if (collisionCheck(pos2, _search_margin)) {
    return false;
  }
  return true;
}

double SkeletonFinder::getNodeRadius(NodePtr curNodePtr) {
  double node_radius = 0;
  for (VertexPtr bv : curNodePtr->black_vertices) {
    node_radius += getDis(curNodePtr->coord, bv->coord);
  }
  for (VertexPtr wv : curNodePtr->white_vertices) {
    node_radius += getDis(curNodePtr->coord, wv->coord);
  }
  node_radius = node_radius / (double)(curNodePtr->black_vertices.size() +
                                       curNodePtr->white_vertices.size());

  return node_radius;
}

double SkeletonFinder::getVerticesRadius(vector<Eigen::Vector3d> vertices) {
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (Eigen::Vector3d v : vertices) {
    center += v;
  }
  center = center / (double)vertices.size();

  double radius = 0;
  for (Eigen::Vector3d v : vertices) {
    radius += getDis(center, v);
  }
  radius = radius / (double)vertices.size();

  return radius;
}

bool SkeletonFinder::checkPtOnFrontier(FrontierPtr ftr_ptr,
                                       Eigen::Vector3d pt) {
  int num_facet = ftr_ptr->facets.size();
  for (int i = 0; i < num_facet; i++) {
    FacetPtr facet = ftr_ptr->facets.at(i);
    if (checkPtOnFacet(facet, pt))
      return true;
    else
      continue;
  }
  return false;
}

bool SkeletonFinder::checkPtOnFacet(FacetPtr facet, Eigen::Vector3d pt) {
  // Take any point on the plane
  Eigen::Vector3d origin = facet->vertices.at(0)->coord;
  Eigen::Vector3d origin_to_pt = pt - origin;
  double dist_signed = origin_to_pt.dot(facet->outwards_unit_normal);

  if (fabs(dist_signed) >= 2 * _search_margin)
    return false;

  Eigen::Vector3d projected_pt = pt - dist_signed * facet->outwards_unit_normal;
  Eigen::Vector3d coord1 = facet->vertices.at(0)->coord;
  Eigen::Vector3d coord2 = facet->vertices.at(1)->coord;
  Eigen::Vector3d coord3 = facet->vertices.at(2)->coord;
  Eigen::Vector3d cross1 = (coord2 - coord1).cross(projected_pt - coord1);
  Eigen::Vector3d cross2 = (coord3 - coord2).cross(projected_pt - coord2);
  Eigen::Vector3d cross3 = (coord1 - coord3).cross(projected_pt - coord3);

  if (cross1(0) * cross2(0) > 0 && cross2(0) * cross3(0) > 0 &&
      cross3(0) * cross1(0) > 0)
    return true;
  else
    return false;
}

FrontierPtr SkeletonFinder::findFlowBackFrontier(Eigen::Vector3d pos,
                                                 int index) {
  for (FrontierPtr f : center_NodeList.at(index)->frontiers) {
    if (getDis(pos, f->proj_center) > 2 * _max_ray_length)
      continue;
    if (checkPtOnFrontier(f, pos))
      return f;
  }

  return NULL;
}

FacetPtr SkeletonFinder::findFlowBackFacet(Eigen::Vector3d pos, int index) {
  for (FacetPtr f : center_NodeList.at(index)->facets) {
    if (getDis(pos, f->center) > 2 * _max_ray_length)
      continue;
    if (checkPtOnFacet(f, pos))
      return f;
  }

  return NULL;
}

bool SkeletonFinder::checkWithinBbx(Eigen::Vector3d pos) {
  return pos(0) >= _x_min && pos(1) >= _y_min && pos(2) >= _z_min &&
         pos(0) <= _x_max && pos(1) <= _y_max && pos(2) <= _z_max;
}

void SkeletonFinder::addBbxToMap(pcl::PointCloud<pcl::PointXYZ> &map) {
  double x_length = _x_max - _x_min;
  double y_length = _y_max - _y_min;
  // double z_length = _z_max - _z_min;
  int x_num = ceil(x_length / _resolution) + 1;
  int y_num = ceil(y_length / _resolution) + 1;
  // int z_num = ceil(z_length / _resolution) + 1;

  if (_is_simulation) { // Add ceiling and floor to search map
    for (int i = 0; i < x_num; i++) {
      for (int j = 0; j < y_num; j++) {
        map.points.push_back(pcl::PointXYZ(_x_min + _resolution * i,
                                           _y_min + _resolution * j, _z_min));
        map.points.push_back(pcl::PointXYZ(_x_min + _resolution * i,
                                           _y_min + _resolution * j, _z_max));
      }
    }
  } else { // Add only ceiling to search map
    for (int i = 0; i < x_num; i++) {
      for (int j = 0; j < y_num; j++) {
        map.points.push_back(pcl::PointXYZ(_x_min + _resolution * i,
                                           _y_min + _resolution * j, _z_max));
      }
    }
  }
}

/* -------------------- visualization ------------------- */

void SkeletonFinder::visualization() {
  visStart();
  visNodesAndVertices();
  visPolygons();
  visFrontiers();
  visMap();
  visConnections();
}

void SkeletonFinder::visNodesAndVertices() {
  // nodes_pcl.clear();
  black_vertices_pcl.clear();
  white_vertices_pcl.clear();
  grey_vertices_pcl.clear();

  int num_nodes = NodeList.size();
  for (int i = num_nodes - 1; i >= 0; i--) {
    if (NodeList.at(i)->rollbacked)
      continue;
    if (NodeList.at(i)->isGate)
      continue;
    // nodes_pcl.points.push_back(pcl::PointXYZ(NodeList.at(i)->coord(0),
    // NodeList.at(i)->coord(1),
    //                                          NodeList.at(i)->coord(2)));

    for (VertexPtr v : NodeList.at(i)->black_vertices) {
      if (v->type == BLACK) {
        black_vertices_pcl.points.push_back(
            pcl::PointXYZ(v->coord(0), v->coord(1), v->coord(2)));
      } else {
        grey_vertices_pcl.points.push_back(
            pcl::PointXYZ(v->coord(0), v->coord(1), v->coord(2)));
      }
    }
    for (VertexPtr v : NodeList.at(i)->white_vertices) {
      white_vertices_pcl.points.push_back(
          pcl::PointXYZ(v->coord(0), v->coord(1), v->coord(2)));
    }

    if (!_visualize_all)
      break;
  }

  pcl::toROSMsg(nodes_pcl, nodes_pcl_ros);
  nodes_pcl_ros.header.frame_id = "map";
  vis_nodes_pub.publish(nodes_pcl_ros);

  pcl::toROSMsg(black_vertices_pcl, black_vertices_pcl_ros);
  black_vertices_pcl_ros.header.frame_id = "map";
  vis_black_vertices_pub.publish(black_vertices_pcl_ros);

  pcl::toROSMsg(white_vertices_pcl, white_vertices_pcl_ros);
  white_vertices_pcl_ros.header.frame_id = "map";
  vis_white_vertices_pub.publish(white_vertices_pcl_ros);

  pcl::toROSMsg(grey_vertices_pcl, grey_vertices_pcl_ros);
  grey_vertices_pcl_ros.header.frame_id = "map";
  vis_grey_vertices_pub.publish(grey_vertices_pcl_ros);
}

void SkeletonFinder::visMap() {
  pcl::toROSMsg(vis_map_pcl, map_pcl_ros);
  map_pcl_ros.header.frame_id = "map";
  map_pcl_ros.height = 1;
  map_pcl_ros.width = vis_map_pcl.points.size();
  vis_map_pub.publish(map_pcl_ros);
}

void SkeletonFinder::visPolygons() {
  visualization_msgs::Marker polygons;

  polygons.header.frame_id = "map";
  polygons.action = visualization_msgs::Marker::ADD;
  polygons.pose.orientation.w = 1.0;
  polygons.type = visualization_msgs::Marker::LINE_LIST;
  polygons.scale.x = 0.03;
  polygons.color.r = 0.25;
  polygons.color.g = 0.75;
  polygons.color.b = 1.0;
  polygons.color.a = 0.2;

  int num_nodes = NodeList.size();

  for (int i = num_nodes - 1; i >= 0; i--) {
    if (NodeList.at(i)->rollbacked)
      continue;
    if (NodeList.at(i)->isGate)
      continue;

    if (_visualize_black_polygon) {
      for (FacetPtr facet : NodeList.at(i)->facets) {
        for (int j = 0; j < 3; j++) {
          geometry_msgs::Point p1, p2;
          p1.x = facet->vertices.at(j)->coord(0);
          p1.y = facet->vertices.at(j)->coord(1);
          p1.z = facet->vertices.at(j)->coord(2);
          p2.x = facet->vertices.at((j + 1) % 3)->coord(0);
          p2.y = facet->vertices.at((j + 1) % 3)->coord(1);
          p2.z = facet->vertices.at((j + 1) % 3)->coord(2);

          polygons.points.push_back(p1);
          polygons.points.push_back(p2);
        }
      }
    } else {
      for (VertexPtr v : NodeList.at(i)->white_vertices) {
        geometry_msgs::Point p1, p2;
        p1.x = v->coord(0);
        p1.y = v->coord(1);
        p1.z = v->coord(2);
        for (VertexPtr v_nbhd : v->connected_vertices) {
          p2.x = v_nbhd->coord(0);
          p2.y = v_nbhd->coord(1);
          p2.z = v_nbhd->coord(2);

          polygons.points.push_back(p1);
          polygons.points.push_back(p2);
        }
      }
      /*
      for (VertexPtr v : NodeList.at(i)->black_vertices) {
        geometry_msgs::Point p1, p2;
        p1.x = v->coord(0);
        p1.y = v->coord(1);
        p1.z = v->coord(2);
        for (VertexPtr v_nbhd : v->connected_vertices) {
          p2.x = v_nbhd->coord(0);
          p2.y = v_nbhd->coord(1);
          p2.z = v_nbhd->coord(2);

          polygons.points.push_back(p1);
          polygons.points.push_back(p2);
        }
      }
      */
    }

    if (!_visualize_all)
      break;
  }

  vis_polygons_pub.publish(polygons);
}

void SkeletonFinder::visFrontiers() {
  visualization_msgs::Marker frontiers;

  frontiers.header.frame_id = "map";
  frontiers.action = visualization_msgs::Marker::ADD;
  frontiers.pose.orientation.w = 1.0;
  frontiers.type = visualization_msgs::Marker::LINE_LIST;
  frontiers.scale.x = 0.03;
  // frontiers.color.r = 0.1;
  // frontiers.color.g = 0.75;
  frontiers.color.g = 1.0;
  // frontiers.color.b = 1.0;
  frontiers.color.a = 0.4;

  int num_nodes = NodeList.size();
  for (int i = num_nodes - 1; i >= 0; i--) {
    if (NodeList.at(i)->rollbacked)
      continue;
    if (NodeList.at(i)->isGate)
      continue;

    for (FrontierPtr ftr : NodeList.at(i)->frontiers) {
      if (ftr->deleted)
        continue;

      // Visualize outwards normal of frontier
      if (_visualize_outwards_normal) {
        geometry_msgs::Point normal;
        normal.x = ftr->outwards_unit_normal(0) + ftr->proj_center(0);
        normal.y = ftr->outwards_unit_normal(1) + ftr->proj_center(1);
        normal.z = ftr->outwards_unit_normal(2) + ftr->proj_center(2);
        geometry_msgs::Point center;
        center.x = ftr->proj_center(0);
        center.y = ftr->proj_center(1);
        center.z = ftr->proj_center(2);

        frontiers.points.push_back(normal);
        frontiers.points.push_back(center);
      }

      for (FacetPtr facet : ftr->facets) {
        /*
        // Visualize facet connections
        if (_visualize_nbhd_facets) {
          int num_nbhd = facet->nbhd_facets.size();
          for (int s = 0; s < num_nbhd; s++) {
            geometry_msgs::Point p1, p2;
            p1.x = facet->center(0);
            p1.y = facet->center(1);
            p1.z = facet->center(2);
            p2.x = facet->nbhd_facets.at(s)->center(0);
            p2.y = facet->nbhd_facets.at(s)->center(1);
            p2.z = facet->nbhd_facets.at(s)->center(2);

            frontiers.points.push_back(p1);
            frontiers.points.push_back(p2);
          }
        }
        */

        // Visualize edges of a facet
        for (int j = 0; j < 3; j++) {
          geometry_msgs::Point p1, p2;
          p1.x = facet->vertices.at(j)->coord(0);
          p1.y = facet->vertices.at(j)->coord(1);
          p1.z = facet->vertices.at(j)->coord(2);
          p2.x = facet->vertices.at((j + 1) % 3)->coord(0);
          p2.y = facet->vertices.at((j + 1) % 3)->coord(1);
          p2.z = facet->vertices.at((j + 1) % 3)->coord(2);

          frontiers.points.push_back(p1);
          frontiers.points.push_back(p2);
        }
      }
    }
    if (!_visualize_all)
      break;
  }
  vis_frontiers_pub.publish(frontiers);
}

void SkeletonFinder::visCurrentFrontier(FrontierPtr ftr) {
  visualization_msgs::Marker frontiers;

  frontiers.header.frame_id = "map";
  // frontiers.action = visualization_msgs::Marker::ADD;
  frontiers.pose.orientation.w = 1.0;
  frontiers.type = visualization_msgs::Marker::LINE_LIST;
  frontiers.scale.x = 0.05;
  frontiers.color.r = 1.0;
  // frontiers.color.g = 1.0;
  // frontiers.color.b = 1.0;
  frontiers.color.a = 0.8;

  geometry_msgs::Point normal;
  normal.x = ftr->outwards_unit_normal(0) + ftr->proj_center(0);
  normal.y = ftr->outwards_unit_normal(1) + ftr->proj_center(1);
  normal.z = ftr->outwards_unit_normal(2) + ftr->proj_center(2);
  geometry_msgs::Point center;
  center.x = ftr->proj_center(0);
  center.y = ftr->proj_center(1);
  center.z = ftr->proj_center(2);

  frontiers.points.push_back(normal);
  frontiers.points.push_back(center);

  for (FacetPtr facet : ftr->facets) {
    for (int j = 0; j < 3; j++) {
      geometry_msgs::Point p1, p2;
      p1.x = facet->vertices.at(j)->coord(0);
      p1.y = facet->vertices.at(j)->coord(1);
      p1.z = facet->vertices.at(j)->coord(2);
      p2.x = facet->vertices.at((j + 1) % 3)->coord(0);
      p2.y = facet->vertices.at((j + 1) % 3)->coord(1);
      p2.z = facet->vertices.at((j + 1) % 3)->coord(2);

      frontiers.points.push_back(p1);
      frontiers.points.push_back(p2);
    }
  }
  vis_cur_frontier_pub.publish(frontiers);
}

void SkeletonFinder::visConnections() {
  visualization_msgs::Marker connections;

  connections.header.frame_id = "map";
  connections.action = visualization_msgs::Marker::ADD;
  connections.pose.orientation.w = 1.0;
  connections.type = visualization_msgs::Marker::LINE_LIST;
  connections.scale.x = 0.25; // 0.075
  connections.color.r = 1.0;
  connections.color.g = 0.392;
  connections.color.b = 0.784;
  connections.color.a = 0.75; // 0.5

  int node_count = 0;
  int connection_count = 0;

  int num_nodes = NodeList.size();
  for (int i = 0; i < num_nodes; i++) {
    NodePtr cur_node = NodeList.at(i);
    if (cur_node->rollbacked)
      continue;
    node_count++;

    if (cur_node->isGate)
      continue;
    int num_connect_nodes = cur_node->connected_Node_ptr.size();
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
      NodePtr connect_node = cur_node->connected_Node_ptr.at(j);

      geometry_msgs::Point p1, p2;
      p1.x = cur_node->coord(0);
      p1.y = cur_node->coord(1);
      p1.z = cur_node->coord(2);
      p2.x = connect_node->coord(0);
      p2.y = connect_node->coord(1);
      p2.z = connect_node->coord(2);

      connections.points.push_back(p1);
      connections.points.push_back(p2);
      connection_count++;
    }
  }

  vis_connections_pub.publish(connections);
  // ROS_WARN("----------------------");
  // ROS_WARN("Eval: Vertices: %d, Edges: %d", node_count, connection_count);
}

void SkeletonFinder::visSphere(Eigen::Vector3d pos, double radius) {
  visualization_msgs::Marker sphere;
  sphere.header.frame_id = "map";
  sphere.action = visualization_msgs::Marker::ADD;
  sphere.type = visualization_msgs::Marker::SPHERE;
  sphere.pose.position.x = pos(0);
  sphere.pose.position.y = pos(1);
  sphere.pose.position.z = pos(2);
  sphere.pose.orientation.w = 1.0;
  sphere.scale.x = 2 * radius;
  sphere.scale.y = 2 * radius;
  sphere.scale.z = 2 * radius;
  sphere.color.r = 1.0;
  sphere.color.b = 1.0;
  sphere.color.a = 0.5;

  vis_spheres_pub.publish(sphere);
}

void SkeletonFinder::visPath() {
  visualization_msgs::Marker vis_path;

  vis_path.header.frame_id = "map";
  vis_path.action = visualization_msgs::Marker::ADD;
  vis_path.pose.orientation.w = 1.0;
  vis_path.type = visualization_msgs::Marker::LINE_LIST;
  vis_path.scale.x = 0.5;
  vis_path.color.r = 0.6;
  vis_path.color.g = 0.196;
  vis_path.color.b = 0.8;
  vis_path.color.a = 1.0;

  int num = path.size();
  for (int i = 0; i < num - 1; i++) {
    Eigen::Vector3d head = path.at(i);
    Eigen::Vector3d tail = path.at(i + 1);

    geometry_msgs::Point p1, p2;
    p1.x = head(0);
    p1.y = head(1);
    p1.z = head(2);
    p2.x = tail(0);
    p2.y = tail(1);
    p2.z = tail(2);

    vis_path.points.push_back(p1);
    vis_path.points.push_back(p2);
  }

  vis_path_pub.publish(vis_path);
}

void SkeletonFinder::visStart() {
  visualization_msgs::Marker marker;

  marker.header.frame_id = "map";
  marker.action = visualization_msgs::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.type = visualization_msgs::Marker::POINTS;
  marker.scale.x = 0.25;
  marker.scale.y = 0.25;
  // marker.color.r = 0.0;
  marker.color.g = 0.0;
  // marker.color.b = 0.0;
  marker.color.a = 1.0;

  geometry_msgs::Point p;
  p.x = _path_start_x;
  p.y = _path_start_y;
  p.z = _path_start_z;

  marker.points.push_back(p);

  p.x = _path_target_x;
  p.y = _path_target_y;
  p.z = _path_target_z;

  marker.points.push_back(p);

  vis_start_pub.publish(marker);
}