#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#include <Eigen/Eigen>
#include <queue>
#include "../utils/quickhull/QuickHull.hpp"
using namespace std;

struct Vertex;
typedef Vertex *VertexPtr;

struct Facet;
typedef Facet *FacetPtr;

struct Frontier;
typedef Frontier *FrontierPtr;

struct Node;
typedef Node *NodePtr;

typedef quickhull::Vector3<double> vec3;
enum vertex_type { WHITE, BLACK, GREY };

struct Vertex {
  Eigen::Vector3d coord;
  Eigen::Vector3d dire_unit_sphere;
  vertex_type type;
  vector<VertexPtr> connected_vertices;
  bool visited;
  // Critical vertices are components of only_black polygon
  bool critical;
  int collision_node_index;
  int sampling_dire_index;
  double dis_to_center;

  Vertex(Eigen::Vector3d coord_, Eigen::Vector3d dire_unit_sphere_, vertex_type type_) {
    coord = coord_;
    dire_unit_sphere = dire_unit_sphere_;
    type = type_;

    // Set default values
    visited = false;
    critical = false;
  }

  Vertex() {}
  ~Vertex() {}
};

struct Facet {
  int index;
  vector<VertexPtr> vertices;

  // Eigen::Vector3d coord1;
  // Eigen::Vector3d coord2;

  Eigen::Vector3d outwards_unit_normal;

  // valid means this facet is contained in some frontiers
  bool valid;
  bool visited;
  bool linked;

  Eigen::Vector3d center;
  Eigen::Vector3d seed_node_pos;

  // Facet plane equation: ax + by + cz + d = 0
  Eigen::Vector4d plane_equation;

  NodePtr master_node;

  vector<FacetPtr> nbhd_facets;
  vector<FacetPtr> nbhd_invalid_facets;

  Facet(vector<VertexPtr> vertices_, NodePtr master_node_) {
    vertices = vertices_;
    master_node = master_node_;

    // Currently, center is the midpoint of the facet
    center = (vertices.at(0)->coord + vertices.at(1)->coord + vertices.at(2)->coord) / 3;

    // Compute plane equation: ax + by + cz + d = 0
    Eigen::Vector3d v1 = vertices.at(1)->coord - vertices.at(0)->coord;
    Eigen::Vector3d v2 = vertices.at(2)->coord - vertices.at(0)->coord;
    Eigen::Vector3d normal = v1.cross(v2);
    normal.normalize();

    double a = normal(0);
    double b = normal(1);
    double c = normal(2);
    Eigen::Vector3d abc;
    abc << a, b, c;
    double d = -abc.dot(vertices.at(0)->coord);

    plane_equation << a, b, c, d;

    // Set default value
    valid = false;
    visited = false;
    linked = false;
  }

  Facet() {}
  ~Facet() {}
};

struct Frontier {
  int index;
  vector<FacetPtr> facets;
  vector<VertexPtr> vertices;

  // Average of the facets center
  Eigen::Vector3d avg_center;
  // Average of the facets outwards_unit_normal;
  Eigen::Vector3d outwards_unit_normal;
  // The projection of avg_coord along outwards_normal,
  // which is on one of the facets owned by the frontier
  Eigen::Vector3d proj_center;
  // Eigen::Vector3d proj_facet_normal;
  // Eigen::Vector3d proj_facet_center;
  FacetPtr proj_facet;

  double cos_theta;

  // Position of the node generate by the frontier
  Eigen::Vector3d next_node_pos;

  // valid means this frontier is clear at the time of processing
  bool valid;
  bool deleted;

  NodePtr master_node;
  NodePtr gate_node;

  Frontier(vector<FacetPtr> facets_, NodePtr master_node_) {
    facets = facets_;
    master_node = master_node_;

    int num_facet = facets.size();
    Eigen::Vector3d coord_sum = Eigen::Vector3d::Zero();
    Eigen::Vector3d normal_sum = Eigen::Vector3d::Zero();
    for (int i = 0; i < num_facet; i++) {
      coord_sum += facets.at(i)->center;
      normal_sum += facets.at(i)->outwards_unit_normal;
    }
    avg_center = coord_sum / num_facet;
    outwards_unit_normal = normal_sum / num_facet;

    // Set default value
    valid = false;
    deleted = false;
    gate_node = NULL;
  }

  Frontier() {}
  ~Frontier() {}
};

struct Node {
  int index;
  Eigen::Vector3d coord;
  Eigen::Vector3d original_coord;
  FrontierPtr seed_frontier;
  double clearance;
  double dis_to_floor;

  bool isGate;
  bool rollbacked;

  vector<VertexPtr> black_vertices;
  vector<VertexPtr> white_vertices;
  // vector<VertexPtr> poly_vertices;
  vector<FacetPtr> facets;
  vector<FrontierPtr> frontiers;
  vector<NodePtr> connected_Node_ptr;

  vector<vec3> sampling_directions;
  vector<vec3> valid_sampling_directions;

  // pcl::search::KdTree<pcl::PointXYZ> polyhedron;

  Node(Eigen::Vector3d coord_, FrontierPtr seed_frontier_, bool isGate_ = false) {
    coord = coord_;
    original_coord = coord_;
    seed_frontier = seed_frontier_;
    isGate = isGate_;

    black_vertices.clear();
    white_vertices.clear();
    facets.clear();
    connected_Node_ptr.clear();

    rollbacked = false;
  }

  Node() {}
  ~Node() {}
};

#endif