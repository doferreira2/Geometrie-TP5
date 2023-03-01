#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Vertex_handle Vertex_handle;

typedef Polyhedron::Point_3 Point;

/// @brief Axis-Aligned bounding box
struct AABB
{
	// TODO...
};

struct OctreeNode
{
	int nb_vertices;
	std::vector<Vertex_handle> vertices;
	std::vector<OctreeNode> children;
};

/// @brief Compute the bounding box of a mesh
/// @param mesh the mesh of interest
/// @return its bounding box
AABB computeBB(const Polyhedron &mesh)
{
	// TODO...
	return AABB{};
}

/// @brief add a level to the given parent Octree node, by creating 8 children with 8 bounding box,
/// sliced in the middle of the parent node
/// @param node the octree node to which 8 children will be added
void addOctreeLevel(OctreeNode &node)
{
	// TODO...
}

constexpr int MAX_POINT = 35; // for testing purposes,
constexpr int MAX_DEPTH = 10; // it would be much better if these values were given to the function where the tree is being constructed.

/// @brief add one vertex to an octree, by following strictly the rules of maximum amount of point in a node, and maximum depth of the tree
/// @param root the root node of the tree
/// @param vert the vertex that will be added, as a Vertex_handle
void addVertexToOctree(OctreeNode &root, Polyhedron::Vertex_handle &vert)
{
	// TODO, this function can be recursive
}

/// @brief A function to generate an octree of the vertices of a mesh,
/// Each vertex will be stored in a node of the octree.
/// the octree shall follow two rules:
///    1- each node shall only contain MAX_POINT vertices
///    2- the depth of tree shall not exceed MAX_DEPTH
///	Remark: the depth of the root is 0
///	Remark: rule 2 wins over rule 1.
/// i.e. a node may contain more vertices than MAX_POINT if the maximum depth is reached.
/// @param mesh the mesh of interest
/// @return an octree node that is the root of the octree for the given mesh

std::vector<Vertex_iterator> getVerticesInBox(Polyhedron &mesh, const Point &boxMin, const Point &boxMax)
{
	std::vector<Vertex_iterator> verticesInBox;

	for (Vertex_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		if (v->point().x() >= boxMin.x() && v->point().y() >= boxMin.y() && v->point().z() >= boxMin.z() && v->point().x() <= boxMax.x() && v->point().y() <= boxMax.y() && v->point().z() <= boxMax.z())
		{
			verticesInBox.push_back(v);
		}
	}

	return verticesInBox;
}

OctreeNode generateOctreeHelper(Polyhedron &mesh, int depth, const Point &min, const Point &max)
{
	OctreeNode node{};

	// if the depth is greater than or equal to the max depth, stop subdividing
	if (depth >= MAX_DEPTH)
	{
		node.vertices = getVerticesInBox(mesh, min, max);
		node.nb_vertices = node.vertices.size();
		return node;
	}

	// if the number of vertices in the box is less than or equal to max point, store the box as a leaf node
	auto vertices = getVerticesInBox(mesh, min, max);
	if (vertices.size() <= MAX_POINT)
	{
		node.vertices = vertices;
		node.nb_vertices = node.vertices.size();
		return node;
	}

	// otherwise, subdivide the box into eight smaller boxes
	Point midpoint = Point(((min.x() + max.x()) / 2.0), ((min.y() + max.y()) / 2.0), ((min.z() + max.z()) / 2.0));

	node.children.resize(8);

	node.children[0] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), min.y(), min.z()), Point(midpoint.x(), midpoint.y(), midpoint.z()));
	node.children[1] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), min.y(), min.z()), Point(max.x(), midpoint.y(), midpoint.z()));
	node.children[2] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), midpoint.y(), min.z()), Point(midpoint.x(), max.y(), midpoint.z()));
	node.children[3] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), midpoint.y(), min.z()), Point(max.x(), max.y(), midpoint.z()));
	node.children[4] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), min.y(), midpoint.z()), Point(midpoint.x(), midpoint.y(), max.z()));
	node.children[5] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), min.y(), midpoint.z()), Point(max.x(), midpoint.y(), max.z()));
	node.children[6] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), midpoint.y(), midpoint.z()), Point(midpoint.x(), max.y(), max.z()));
	node.children[7] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), midpoint.y(), midpoint.z()), Point(max.x(), max.y(), max.z()));

	return node;
}

OctreeNode generateOctree(Polyhedron &mesh /*, max number of point, max depth...*/)
{
	// start by defining the bounding box of the mesh
	Point min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	Point max(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

	for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		min = Point(std::min(min.x(), v->point().x()), std::min(min.y(), v->point().y()), std::min(min.z(), v->point().z()));
		max = Point(std::max(max.x(), v->point().x()), std::max(max.y(), v->point().y()), std::max(max.z(), v->point().z()));
	}

	return generateOctreeHelper(mesh, 0, min, max);
}

/// @brief find a specific vertex inside an octree (using a dichotomy algorithm)
/// @param vh the vertex handle to look for
/// @return the address of the node (not the prettiest way, feel free to handle it differently)
OctreeNode *findVertexInOctree(OctreeNode &root, Polyhedron::Vertex_const_handle &vh)
{
	// TODO....
	return &root;
}

/// @brief (optional) Utility function that takes an octree and apply a function (or more useful, a lambda !)
/// to each leaf of the Octree (each node containing vertices).
/// Can be useful to avoid declaring a new recursive function each time...
/// @param root the root node of the Octree of interest
/// @param func a lambda supposed to do something on a given Octree node.
void browseNodes(const OctreeNode &root, std::function<void(const OctreeNode &)> func)
{
	// if there are no vertices in the node we do nothing

	// if the nodes contains vertices, then "func" is called on the node

	// go through all the children of the current node and calls browseNodes recursively.
	// browseNodes(/*TODO*/, func);
}

void extractMeshFromOctree(const OctreeNode &root, const Polyhedron &mesh)
{

	std::vector<Point> vertices;
	std::vector<std::vector<int>> faces;

	// TODO: fill "vertices" and "faces" by going through the octree

	std::ofstream out("octree_meshres.off");
	out << "OFF" << std::endl;
	out << vertices.size() << " " << faces.size() << " 0" << std::endl;
	for (const auto &v : vertices)
	{
		out << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	for (const auto &f : faces)
	{
		out << f.size() << " ";
		for (auto fi : f)
		{
			out << fi << " ";
		}
		out << std::endl;
	}
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}

	unsigned int nbVerts = 0;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		++nbVerts;
	}
	std::cout << "Nombre de sommets: " << nbVerts << std::endl;

	unsigned int nbEdges = 0;
	for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		++nbEdges;
	}
	nbEdges /= 2;
	std::cout << "Nombre d'arêtes: " << nbEdges << std::endl;

	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		++nbFaces;
	}
	std::cout << "Nombre de faces: " << nbFaces << std::endl;

	const auto octree = generateOctree(mesh);

	extractMeshFromOctree(octree, mesh);

	return 0;
}
