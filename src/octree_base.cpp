#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_const_iterator Facet_iterator;
typedef Polyhedron::Vertex_const_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_const_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator Halfedge_facet_circulator;

// typedef Polyhedron::Vertex_handle Vertex_handle;

typedef Polyhedron::Point_3 Point;

typedef struct s_color
{
	float R;
	float V;
	float B;
} color;

color colorPalette[8] = {
	{0, 0, 0},
	{0, 0, 255},
	{0, 255, 0},
	{0, 255, 255},
	{255, 0, 0},
	{255, 0, 255},
	{255, 255, 0},
	{255, 255, 255},
};

/// @brief Axis-Aligned bounding box
struct AABB
{
	Point boxMin;
	Point boxMax;

	// TODO...
};

struct OctreeNode
{
	int nb_vertices;
	int profondeur;
	int nieme;
	AABB bbox;
	std::vector<Vertex_iterator> vertices;
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
constexpr int MAX_DEPTH = 5;  // it would be much better if these values were given to the function where the tree is being constructed.

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

std::vector<Vertex_iterator> getVerticesInBox(const Polyhedron &mesh, const Point &boxMin, const Point &boxMax)
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

OctreeNode generateOctreeHelper(const Polyhedron &mesh, int depth, const Point &min, const Point &max, int i)
{

	OctreeNode node{};
	node.profondeur = depth;
	node.bbox.boxMax = max;
	node.bbox.boxMin = min;
	node.nieme = i;

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

	node.children[0] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), min.y(), min.z()), Point(midpoint.x(), midpoint.y(), midpoint.z()), 0);
	node.children[1] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), min.y(), min.z()), Point(max.x(), midpoint.y(), midpoint.z()), 1);
	node.children[2] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), midpoint.y(), min.z()), Point(midpoint.x(), max.y(), midpoint.z()), 2);
	node.children[3] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), midpoint.y(), min.z()), Point(max.x(), max.y(), midpoint.z()), 3);
	node.children[4] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), min.y(), midpoint.z()), Point(midpoint.x(), midpoint.y(), max.z()), 4);
	node.children[5] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), min.y(), midpoint.z()), Point(max.x(), midpoint.y(), max.z()), 5);
	node.children[6] = generateOctreeHelper(mesh, depth + 1, Point(min.x(), midpoint.y(), midpoint.z()), Point(midpoint.x(), max.y(), max.z()), 6);
	node.children[7] = generateOctreeHelper(mesh, depth + 1, Point(midpoint.x(), midpoint.y(), midpoint.z()), Point(max.x(), max.y(), max.z()), 7);

	return node;
}

OctreeNode generateOctree(const Polyhedron &mesh /*, max number of point, max depth...*/)
{
	// start by defining the bounding box of the mesh
	Point max(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	Point min(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

	for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		min = Point(std::min(min.x(), v->point().x()), std::min(min.y(), v->point().y()), std::min(min.z(), v->point().z()));
		max = Point(std::max(max.x(), v->point().x()), std::max(max.y(), v->point().y()), std::max(max.z(), v->point().z()));
	}

	return generateOctreeHelper(mesh, 0, min, max, 0);
}

bool VertexInBBbox(AABB &box, Vertex_iterator &v)
{

	if (v->point().x() >= box.boxMin.x() && v->point().y() >= box.boxMin.y() && v->point().z() >= box.boxMin.z() && v->point().x() <= box.boxMax.x() && v->point().y() <= box.boxMax.y() && v->point().z() <= box.boxMax.z())
	{

		return true;
	}
	else
	{
		return false;
	}
}

/// @brief find a specific vertex inside an octree (using a dichotomy algorithm)
/// @param vh the vertex handle to look for
/// @return the address of the node (not the prettiest way, feel free to handle it differently)
OctreeNode *findVertexInOctree(OctreeNode &root, Vertex_iterator &vh)
{

	// Si le nœud courant ne contient pas la boîte englobante du sommet, il n'y a pas de point
	// correspondant dans cet octree, on retourne donc nullptr
	if (!VertexInBBbox(root.bbox, vh))
	{
		return nullptr;
	}

	// Si le nœud courant contient le sommet, on le retourne
	if (std::find(root.vertices.begin(), root.vertices.end(), vh) != root.vertices.end()) // root.vertices. bo find(vh) != root.points.end())
	{
		return &root;
	}

	// Si le nœud courant ne contient pas le sommet, on recherche récursivement
	// dans les nœuds enfants qui contiennent la boîte englobante du sommet
	for (OctreeNode &child : root.children)
	{
		if (VertexInBBbox(child.bbox, vh))
		{
			OctreeNode *result = findVertexInOctree(child, vh);
			if (result != nullptr)
			{
				return result;
			}
		}
	}

	// Si le sommet n'a pas été trouvé dans l'octree, on retourne nullptr
	return nullptr;
}

/// @brief (optional) Utility function that takes an octree and apply a function (or more useful, a lambda !)
/// to each leaf of the Octree (each node containing vertices).
/// Can be useful to avoid declaring a new recursive function each time...
/// @param root the root node of the Octree of interest
/// @param func a lambda supposed to do something on a given Octree node.
void browseNodes(OctreeNode &root, std::function<void(const OctreeNode &)> func)
{
	if (root.nb_vertices != 0)
	{
		func(root);
		return;
	}

	for (OctreeNode &child : root.children)
	{
		browseNodes(child, func);
	}

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

void writeJSONfromOctree(const OctreeNode &tree, std::ofstream &file)
{
	int i = tree.nieme;
	file << " { \"Profondeur \" : " << tree.profondeur << "," << std::endl;
	file << "  \"vertex \" : " << tree.nb_vertices << "," << std::endl;
	file << "  \"i \" : " << i << "," << std::endl;
	file << "  \"enfant \" : [ " << std::endl;

	if (tree.profondeur == 0)
		i = 7;

	for (auto &t : tree.children)
	{

		writeJSONfromOctree(t, file);
		// i = (i >= 7 ? 0 : i + 1);
	}

	file << "]" << std::endl;
	file << (i >= 7 ? "}" : "},") << std::endl;
	//file << "  \"i \" : " << i << "," << std::endl;
}

void writeCOFFfromMeshOctree(const Polyhedron &mesh, OctreeNode &tree, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
			  << mesh.size_of_vertices() << ' '
			  << mesh.size_of_facets() << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

	for (Vertex_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		in_myfile << v->point();

		auto *noeud = findVertexInOctree(tree, v);

		if (noeud != nullptr)
		{
			std::cout << "prof : " << noeud->profondeur << std::endl;

			auto redValue = colorPalette[noeud->profondeur].R;
			auto greenValue = colorPalette[noeud->profondeur].V;
			auto blueValue = colorPalette[noeud->profondeur].B;

			in_myfile << " " << redValue << " " << greenValue << " " << blueValue;
		}
		in_myfile << std::endl;
	}

	// std::copy(mesh.points_begin(), mesh.points_end(),std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::setprecision(5) << std::fixed; // set the format of floats to X.XXXXX

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

void ComputPointOctree(std::vector<Point> &vect, OctreeNode &tree, std::ofstream &file)
{

	if (tree.nb_vertices != 0)
	{
		vect.push_back(Point(tree.bbox.boxMax.x(), tree.bbox.boxMax.y(), tree.bbox.boxMax.z()));
		vect.push_back(Point(tree.bbox.boxMax.x(), tree.bbox.boxMax.y(), tree.bbox.boxMin.z()));
		vect.push_back(Point(tree.bbox.boxMax.x(), tree.bbox.boxMin.y(), tree.bbox.boxMax.z()));
		vect.push_back(Point(tree.bbox.boxMax.x(), tree.bbox.boxMin.y(), tree.bbox.boxMin.z()));
		vect.push_back(Point(tree.bbox.boxMin.x(), tree.bbox.boxMax.y(), tree.bbox.boxMax.z()));
		vect.push_back(Point(tree.bbox.boxMin.x(), tree.bbox.boxMax.y(), tree.bbox.boxMin.z()));
		vect.push_back(Point(tree.bbox.boxMin.x(), tree.bbox.boxMin.y(), tree.bbox.boxMax.z()));
		vect.push_back(Point(tree.bbox.boxMin.x(), tree.bbox.boxMin.y(), tree.bbox.boxMin.z()));
	}
	else
	{
		for (OctreeNode &child : tree.children)
		{
			ComputPointOctree(vect, child, file);
		}
	}
}

void writeFaceHelper(int nbcube, std::ofstream &file)
{
	for (int i = 0; i < nbcube; i++)
	{
		file << '4' << ' ' << 0 + (8 * i) << ' ' << 4 + (8 * i) << ' ' << 6 + (8 * i) << ' ' << 2 + (8 * i) << std::endl;
		file << '4' << ' ' << 3 + (8 * i) << ' ' << 2 + (8 * i) << ' ' << 6 + (8 * i) << ' ' << 7 + (8 * i) << std::endl;
		file << '4' << ' ' << 7 + (8 * i) << ' ' << 6 + (8 * i) << ' ' << 4 + (8 * i) << ' ' << 5 + (8 * i) << std::endl;
		file << '4' << ' ' << 5 + (8 * i) << ' ' << 1 + (8 * i) << ' ' << 3 + (8 * i) << ' ' << 7 + (8 * i) << std::endl;
		file << '4' << ' ' << 1 + (8 * i) << ' ' << 0 + (8 * i) << ' ' << 2 + (8 * i) << ' ' << 3 + (8 * i) << std::endl;
		file << '4' << ' ' << 5 + (8 * i) << ' ' << 4 + (8 * i) << ' ' << 0 + (8 * i) << ' ' << 1 + (8 * i) << std::endl;
	}
}

void writeOFFfromOctree(OctreeNode &tree, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	std::vector<Point> v_point;

	std::cout << "init : " << v_point.size() << std::endl;

	ComputPointOctree(v_point, tree, in_myfile);

	std::cout << v_point.size() << std::endl;

	in_myfile << "OFF" << std::endl // "COFF" makes the file support color informations
			  << v_point.size() << ' '
			  << (6 * (v_point.size() / 8)) << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(v_point.begin(), v_point.end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	writeFaceHelper((v_point.size() / 8), in_myfile);

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

// TC : jr sjui dorian, je mange des cartes arduinis au petit dej, maim miam les pcb vive l'eltricté, paul pinault le boss, je veux lui faire des choses

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

	auto octree = generateOctree(mesh);

	// extractMeshFromOctree(octree, mesh);

	writeOFFfromOctree(octree, "Octree.off");

	std::ofstream file;
	file.open("Octree.json");
	file << "{" << std::endl;
	file << "\"meshName\": \" " << argv[1] << "\",\"tree\": [ " << std::endl;
	writeJSONfromOctree(octree, file);
	file << "]}" << std::endl;

	file.close();

	return 0;
}
