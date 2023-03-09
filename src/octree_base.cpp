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

typedef Polyhedron::Point_3 Point;
typedef std::vector<Vertex_iterator> Vertex_vector;

typedef struct s_color
{
	float R;
	float V;
	float B;
} color;

color colorPalette[8] = {
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 255},
	{0.0, 255, 0.0},
	{0.0, 255, 255},
	{255, 0.0, 0.0},
	{255, 0.0, 255},
	{255, 255, 0.0},
	{255, 255, 255},
};

/// @brief Axis-Aligned bounding box
struct AABB
{
	Point boxMin;
	Point boxMax;
};

struct OctreeNode
{
	int nb_vertices;
	int profondeur;
	int nieme;
	AABB bbox;
	Vertex_vector vertices;
	std::vector<OctreeNode> children;
};

int MAX_POINT; // for testing purposes,
int MAX_DEPTH; // it would be much better if these values were given to the function where the tree is being constructed.

typedef std::map<OctreeNode *, int> node_int_map;

Vertex_vector getVerticesInBox(const Polyhedron &mesh, const Point &boxMin, const Point &boxMax)
{
	Vertex_vector verticesInBox;

	for (Vertex_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		if (v->point().x() >= boxMin.x() && v->point().y() >= boxMin.y() && v->point().z() >= boxMin.z() && v->point().x() <= boxMax.x() && v->point().y() <= boxMax.y() && v->point().z() <= boxMax.z())
		{
			verticesInBox.push_back(v);
		}
	}

	return verticesInBox;
}

OctreeNode generateOctreeHelper(const Polyhedron &mesh, int depth, const Point &min, const Point &max, const int i)
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
	Vertex_vector vertices = getVerticesInBox(mesh, min, max);
	if (vertices.size() <= MAX_POINT)
	{
		node.vertices = vertices;
		node.nb_vertices = node.vertices.size();
		return node;
	}

	// otherwise, subdivide the box into eight smaller boxes
	const Point midpoint = Point(((min.x() + max.x()) / 2.0), ((min.y() + max.y()) / 2.0), ((min.z() + max.z()) / 2.0));

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
OctreeNode generateOctree(const Polyhedron &mesh)
{
	// start by defining the bounding box of the mesh
	Point min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	Point max(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

	for (Vertex_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
	{
		min = Point(std::min(min.x(), v->point().x()), std::min(min.y(), v->point().y()), std::min(min.z(), v->point().z()));
		max = Point(std::max(max.x(), v->point().x()), std::max(max.y(), v->point().y()), std::max(max.z(), v->point().z()));
	}

	return generateOctreeHelper(mesh, 0, min, max, 0);
}

bool VerticeInBox(const AABB &box, const Vertex_iterator &v)
{
	if (v->point().x() >= box.boxMin.x() && v->point().y() >= box.boxMin.y() && v->point().z() >= box.boxMin.z() && v->point().x() <= box.boxMax.x() && v->point().y() <= box.boxMax.y() && v->point().z() <= box.boxMax.z())
	{
		return true;
	}
	return false;
}

/// @brief find a specific vertex inside an octree (using a dichotomy algorithm)
/// @param vh the vertex handle to look for
/// @return the address of the node (not the prettiest way, feel free to handle it differently)
OctreeNode *findVertexInOctree(OctreeNode &root, const Vertex_iterator &vh)
{

	// Si le nœud courant ne contient pas la boîte englobante du sommet, il n'y a pas de point
	// correspondant dans cet octree, on retourne donc nullptr
	if (!VerticeInBox(root.bbox, vh))
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
		if (VerticeInBox(child.bbox, vh))
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
	// if the nodes contains vertices, then "func" is called on the node
	if (root.nb_vertices != 0)
	{
		func(root);
		return;
	}
	// go through all the children of the current node and calls browseNodes recursively.

	for (OctreeNode &child : root.children)
	{
		browseNodes(child, func);
	}

	// if there are no vertices in the node we do nothing
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

	for (const OctreeNode &t : tree.children)
	{
		writeJSONfromOctree(t, file);
		// i = (i >= 7 ? 0 : i + 1);
	}

	file << "]" << std::endl;
	file << (i >= 7 ? "}" : "},") << std::endl;
	// file << "  \"i \" : " << i << "," << std::endl;
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

		OctreeNode *noeud = findVertexInOctree(tree, v);

		if (noeud != nullptr)
		{
			float redValue = colorPalette[noeud->nieme].R / ((MAX_DEPTH + 1) - noeud->profondeur);
			float greenValue = colorPalette[noeud->nieme].V / ((MAX_DEPTH + 1) - noeud->profondeur);
			float blueValue = colorPalette[noeud->nieme].B / ((MAX_DEPTH + 1) - noeud->profondeur);

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

void ComputPointOctree(std::vector<Point> &vect, const OctreeNode &tree)
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
		for (const OctreeNode &child : tree.children)
		{
			ComputPointOctree(vect, child);
		}
	}
}

void writeFaceHelper(const int nbcube, std::ofstream &file)
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

void writeOFFfromOctree(const OctreeNode &tree, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	std::vector<Point> v_point;

	ComputPointOctree(v_point, tree);

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

void computSimplePoint(std::vector<Point> &vect, node_int_map &ind, OctreeNode &tree)
{
	if (tree.nb_vertices != 0)
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for (const Vertex_iterator &it : tree.vertices)
		{
			x += it->point().x();
			y += it->point().y();
			z += it->point().z();
		}
		vect.push_back(Point((x / tree.nb_vertices), (y / tree.nb_vertices), (z / tree.nb_vertices)));
		ind[&tree] = (vect.size() - 1);
	}
	else
	{
		for (OctreeNode &child : tree.children)
		{
			computSimplePoint(vect, ind, child);
		}
	}
}

void simplifMesh(OctreeNode &tree, Polyhedron &mesh, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	std::ostringstream fil;

	CGAL::set_ascii_mode(in_myfile);

	std::vector<Point> v_point;
	node_int_map iM_node;

	computSimplePoint(v_point, iM_node, tree);
	int nb_face = 0;
	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		// get the three vertices of the current face
		Vertex_iterator v1 = f->halfedge()->vertex();
		OctreeNode *node1 = findVertexInOctree(tree, v1);
		// get the nodes of the octree containing each vertex

		Vertex_iterator v2 = f->halfedge()->next()->vertex();
		OctreeNode *node2 = findVertexInOctree(tree, v2);

		Vertex_iterator v3 = f->halfedge()->next()->next()->vertex();
		OctreeNode *node3 = findVertexInOctree(tree, v3);

		if (f->is_quad())
		{
			Vertex_iterator v4 = f->halfedge()->next()->next()->next()->vertex();
			OctreeNode *node4 = findVertexInOctree(tree, v4);

			if (node1 != node2 && node1 != node3 && node1 != node4 && node2 != node3 && node2 != node4 && node3 != node4)
			{
				fil << "4 " << iM_node[node1] << " " << iM_node[node2] << " " << iM_node[node3] << " " << iM_node[node4] << " \n";
				nb_face++;
			}
		}
		else
		{
			if (node1 != node2 && node1 != node3 && node2 != node3)
			{
				fil << "3 " << iM_node[node1] << " " << iM_node[node2] << " " << iM_node[node3] << " \n";
				nb_face++;
			}
		}

		// check if all three vertices are in diffrent node
	}

	in_myfile << "OFF" << std::endl // "COFF" makes the file support color informations
			  << v_point.size() << ' '
			  << nb_face << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(v_point.begin(), v_point.end(), std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	in_myfile << fil.str();

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

	MAX_DEPTH = (argc >= 3 ? atoi(argv[2]) : 10);
	MAX_POINT = (argc >= 4 ? atoi(argv[3]) : 50);

	OctreeNode octree = generateOctree(mesh);

	std::ofstream file;
	file.open("Octree.json");
	file << "{" << std::endl;
	file << "\"meshName\": \" " << argv[1] << "\",\"tree\": [ " << std::endl;
	writeJSONfromOctree(octree, file);
	file << "]}" << std::endl;
	file.close();

	writeCOFFfromMeshOctree(mesh, octree, "colorMesh.off");
	writeOFFfromOctree(octree, "Octree.off");

	simplifMesh(octree, mesh, "simple.off");

	// simplifi(mesh, octree);
	return 0;
}

/* BROUILLON

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


*/