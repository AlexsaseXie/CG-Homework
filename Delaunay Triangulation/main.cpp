#include <cstdio>
#include <algorithm>
using namespace std;

const int POINT_MAX = 1e7 + 5;
typedef long long ll;
const double eps = 0;

class DAGNode;
class Point2D;
class Edge;
class Triangle;

class Point2D {
public:
	int x;
	int y;
	int true_id;
public:
	Point2D() { x = 0; y = 0; }
	Point2D(int mx, int my) { x = mx; y = my; }
	~Point2D() {}
};

typedef Point2D Point;

class Edge {
public:
	int start_vertex_index;
	Edge * twin;
	Triangle * inc;
public:
	bool counted = false;
public:
	Edge() { start_vertex_index = -1; twin = nullptr; inc = nullptr; }
	Edge(int a) { start_vertex_index = a; twin = nullptr; inc = nullptr; }
	~Edge() {}
public:
	void set_twin(Edge * tw) { twin = tw; }
	void set_inc(Triangle * in) { inc = in; }
	void set_start_vertex_index(int a) { start_vertex_index = a; }
};

void determine_center_radius(Point &p1, Point &p2, Point &p3, double &center_x, double &center_y, double &square_radius) {
	ll A1 = 2 * (p2.x - p1.x);
	ll B1 = 2 * (p2.y - p1.y);
	double C1 = (ll)p2.x * p2.x + (ll)p2.y * p2.y - (ll)p1.x * p1.x - (ll)p1.y * p1.y;
	ll A2 = 2 * (p3.x - p2.x);
	ll B2 = 2 * (p3.y - p2.y);
	double C2 = (ll)p3.x * p3.x + (ll)p3.y * p3.y - (ll)p2.x * p2.x - (ll)p2.y * p2.y;
	double K1 = (C1 * (double)B2) - (C2 * (double)B1);
	double K2 = ((double)A1 * C2) - ((double)A2 * C1);
	ll DIV = (A1 * B2) - (A2 * B1);

	center_x = (double)K1 / (double)DIV;
	center_y = (double)K2 / (double)DIV;
	square_radius = ((double)p1.x - center_x) * ((double)p1.x - center_x) + ((double)p1.y - center_y) * ((double)p1.y - center_y);
}

class DAGNode {
public:
	Triangle * tri;
	DAGNode * children[3];
	bool visited = false;
public:
	DAGNode() { tri = nullptr; children[0] = nullptr; children[1] = nullptr; children[2] = nullptr;}
	~DAGNode() {}
public:
	bool is_leaf() {
		return (children[0] == nullptr && children[1] == nullptr && children[2] == nullptr);
	}
};

class Triangle {
public:
	int ids[3];	// counter-clockwise
	Edge edges[3]; // correspond to ids; containing adja info
	DAGNode * dag_node = nullptr;
	double center_x;
	double center_y;
	double square_radius = -1;
public:
	void init(int a, int b, int c) {
		ids[0] = a; ids[1] = b; ids[2] = c;
		for (int i = 0; i < 3; i++) {
			edges[i].set_start_vertex_index(i);
			edges[i].set_inc(this);
		}
		dag_node = nullptr;
	}
	Triangle() { 
		ids[0] = -1; ids[1] = -1; ids[2] = -1;
		dag_node = nullptr;
	}
	Triangle(int a, int b, int c) { 
		init(a, b, c);
	}
	Triangle(int a, int b, int c, Point * pts) {
		init(a, b, c);
		//determine_center_radius(pts[ids[0]], pts[ids[1]], pts[ids[2]], center_x, center_y, square_radius);
	}
	~Triangle() { }
};

ll cross_product(int p1x, int p1y, int p2x, int p2y) {
	return ((long long)p1x * p2y) - ((long long)p1y * p2x);
}

ll to_left_value(Point a, Point b, Point c) {
	int b_a_x = b.x - a.x;
	int b_a_y = b.y - a.y;
	int c_a_x = c.x - a.x;
	int c_a_y = c.y - a.y;
	ll cross_pro = cross_product(b_a_x, b_a_y, c_a_x, c_a_y);
	return cross_pro;
}

bool to_left(Point a, Point b, Point c) {
	ll cross_pro = to_left_value(a, b, c);
	return cross_pro >= 0;
}

bool inside_circle(Point &p1, Triangle * current_tri, Point * pts) {
	if (current_tri->square_radius == -1) {
		determine_center_radius(pts[current_tri->ids[0]], pts[current_tri->ids[1]], pts[current_tri->ids[2]], current_tri->center_x, current_tri->center_y, current_tri->square_radius);
	}

	double p_square_radius = ((double)p1.x - current_tri->center_x) * ((double)p1.x - current_tri->center_x)
		+ ((double)p1.y - current_tri->center_y) * ((double)p1.y - current_tri->center_y);
	return p_square_radius - current_tri->square_radius < - eps;
}

bool point_greater(Point &a, Point &b) {
	if (a.y > b.y) return true;
	if (a.y < b.y) return false;
	return a.x > b.x;
}

bool inside_triangle(Point &pr, Triangle *tri, Point * pts, int n) {
	ll to_left_values[3];
	for (int i = 0; i < 3; i++) {
		int id1 = tri->ids[i];
		int id2 = tri->ids[(i + 1) % 3];
		if (id1 >= n && id2 >= n) {  // p-2 p-1
			to_left_values[i] = 1; 
			continue; 
		}
		if (id1 >= n) { // p-1 pi or p-2 pi
			if (id1 == n) {
				// p-1 
				bool flag = point_greater(pr, pts[id2]);
				if (!flag) to_left_values[i] = 1;
				else to_left_values[i] = -1;
			} 
			else {
				// p-2 
				bool flag = point_greater(pr, pts[id2]);
				if (flag) to_left_values[i] = 1;
				else to_left_values[i] = -1;
			}
			continue; 
		}
		if (id2 >= n) { // pi p-1 or pi p-2
			if (id2 == n) {
				// p-1 
				bool flag = point_greater(pr, pts[id1]);
				if (flag) to_left_values[i] = 1;
				else to_left_values[i] = -1;
			}
			else {
				// p-2 
				bool flag = point_greater(pr, pts[id1]);
				if (!flag) to_left_values[i] = 1;
				else to_left_values[i] = -1;
			}
			continue;
		}
		to_left_values[i] = to_left_value(pts[id1], pts[id2], pr);
	}

	return (to_left_values[0] >= 0 && to_left_values[1] >= 0 && to_left_values[2] >= 0);
}

DAGNode * find_triangle_contain_pr(DAGNode * head, Point &pr, Point * pts, int n) {
	if (head->is_leaf())
		return head;

	for (int i = 0; i < 3; i++) {
		if (inside_triangle(pr, head->children[i]->tri, pts, n)) {
			return find_triangle_contain_pr(head->children[i], pr, pts, n);
		}
	}

	// failed
	return nullptr;
}

bool is_convex_quad(int pr, int pk, int pi, int pj, Point * pts, int n) {
	// judge whether pr pi pk pj is convex quad

	if (pk < n) {
		ll to_left_value_i, to_left_value_j;
		int prk_low, prk_high;
		if (point_greater(pts[pr], pts[pk])) {
			prk_low = pk;
			prk_high = pr;
		}
		else {
			prk_low = pr;
			prk_high = pk;
		}

		if (pi >= n && pj >= n)
			return false;

		if (pi == n) to_left_value_i = -1; 
		else {
			if (pi == n + 1) to_left_value_i = 1;
			else 
				to_left_value_i = to_left_value(pts[prk_low], pts[prk_high], pts[pi]);
		}

		if (pj == n) to_left_value_j = -1;
		else {
			if (pj == n + 1) to_left_value_j = 1;
			else 
				to_left_value_j = to_left_value(pts[prk_low], pts[prk_high], pts[pj]);
		}

		return (to_left_value_i < 0 && to_left_value_j > 0) || (to_left_value_i > 0 && to_left_value_j < 0);
	}
	else {
		// pk >= n
		if (pi >= n || pj >= n) return false;
		else
			return false;
	}
}

void legalize(int pr_in_tri_index, int pi_in_tri_index, int pj_in_tri_index, DAGNode * pi_pj_pr_dagnode, Point * pts, int n) {
	Triangle * pi_pj_pr_tri = pi_pj_pr_dagnode->tri;
	Edge* pi_pj = &pi_pj_pr_tri->edges[pi_in_tri_index];
	Edge* adja_edge_pj_pi = pi_pj->twin;

	if (adja_edge_pj_pi == nullptr) return;

	int pj_pi_pk_tri_pj_index = adja_edge_pj_pi->start_vertex_index;
	int pj_pi_pk_tri_pi_index = (pj_pi_pk_tri_pj_index + 1) % 3;
	int pj_pi_pk_tri_pk_index = (pj_pi_pk_tri_pj_index + 2) % 3;
	Triangle * pj_pi_pk_tri = adja_edge_pj_pi->inc;
	DAGNode * pj_pi_pk_dagnode = pj_pi_pk_tri->dag_node;

	int pk = pj_pi_pk_tri->ids[pj_pi_pk_tri_pk_index];
	int pr = pi_pj_pr_tri->ids[pr_in_tri_index];
	int pi = pi_pj_pr_tri->ids[pi_in_tri_index];
	int pj = pi_pj_pr_tri->ids[pj_in_tri_index];

	bool notlegal = true;
	// p0 p-1 p-2 legal
	if (pi >= n - 1 && pj >= n - 1) {
		notlegal = false;
		return;
	}

	if (pi >= n || pj >= n || pk >= n) {
		int tmp_pk = pk >= n - 1 ? -(pk - (n - 1)) : pk + 1;
		int tmp_pr = pr >= n - 1 ? -(pr - (n - 1)) : pr + 1;
		int tmp_pi = pi >= n - 1 ? -(pi - (n - 1)) : pi + 1;
		int tmp_pj = pj >= n - 1 ? -(pj - (n - 1)) : pj + 1;
		int min_r_k = tmp_pk < tmp_pr ? tmp_pk : tmp_pr;
		int min_i_j = tmp_pi < tmp_pj ? tmp_pi : tmp_pj;
		notlegal = !(min_r_k < min_i_j);
	}
	else {
		notlegal = inside_circle(pts[pk], pi_pj_pr_tri, pts);
	}
	if (notlegal) {
		if (!is_convex_quad(pr, pk, pi, pj, pts, n))
			return;

		// flip
		int ikj_vetices[3]= { pi_pj_pr_tri->ids[pi_in_tri_index], pk , pi_pj_pr_tri->ids[pj_in_tri_index] };
		
		DAGNode * new_dagnode[2];
		Triangle * new_tris[2];

		for (int i = 0; i < 2; i++) {
			new_dagnode[i] = new DAGNode;
			new_tris[i] = new Triangle(ikj_vetices[i], ikj_vetices[i + 1], pr, pts);

			new_dagnode[i]->tri = new_tris[i];
			new_tris[i]->dag_node = new_dagnode[i];
		}

		//outer_adja_edge[4]
		Edge * outer_adja_edge[4] = {
			pi_pj_pr_tri->edges[pr_in_tri_index].twin,
			pj_pi_pk_tri->edges[pj_pi_pk_tri_pi_index].twin,
			pj_pi_pk_tri->edges[pj_pi_pk_tri_pk_index].twin,
			pi_pj_pr_tri->edges[pj_in_tri_index].twin 
		};
		Edge * new_tris_outer_edge[4] = {
			&new_tris[0]->edges[2],
			&new_tris[0]->edges[0],
			&new_tris[1]->edges[0],
			&new_tris[1]->edges[1]
		};
		for (int i = 0; i < 4; i++) {
			if (outer_adja_edge[i] != nullptr) {
				outer_adja_edge[i]->twin = new_tris_outer_edge[i];
			}
			new_tris_outer_edge[i]->twin = outer_adja_edge[i];
		}

		new_tris[0]->edges[1].twin = &new_tris[1]->edges[2];
		new_tris[1]->edges[2].twin = &new_tris[0]->edges[1];

		pi_pj_pr_dagnode->children[0] = new_dagnode[0];
		pi_pj_pr_dagnode->children[1] = new_dagnode[1];
		pj_pi_pk_dagnode->children[0] = new_dagnode[0];
		pj_pi_pk_dagnode->children[1] = new_dagnode[1];


		legalize(2, 0, 1, new_dagnode[0], pts, n);
		legalize(2, 0, 1, new_dagnode[1], pts, n);
	}

	return;
}

void insert_pr_inside(int pr, DAGNode * pi_pj_pk_dagnode, Point * pts) {
	// insert 3 triangles
	Triangle * pi_pj_pk_tri = pi_pj_pk_dagnode->tri;
	Triangle * new_tris[3];
	for (int i = 0; i < 3; i++) {
		// new dagnode & tris
		pi_pj_pk_dagnode->children[i] = new DAGNode;
		new_tris[i] = new Triangle(pi_pj_pk_tri->ids[i], pi_pj_pk_tri->ids[(i + 1) % 3], pr, pts);

		pi_pj_pk_dagnode->children[i]->tri = new_tris[i];
		new_tris[i]->dag_node = pi_pj_pk_dagnode->children[i];

		Edge * pi_pj_pk_current_edge = & pi_pj_pk_tri->edges[i];
		// find adja tri
		Edge * adja_edge = pi_pj_pk_current_edge->twin;
		
		// link outer adja pointer
		if (adja_edge != nullptr) {
			adja_edge->twin = &new_tris[i]->edges[0];
		}
		new_tris[i]->edges[0].twin = adja_edge;
	}

	// link inner adja pointer
	for (int i = 0; i < 3; i++) {
		new_tris[i]->edges[1].twin = &(new_tris[(i + 1) % 3]->edges[2]);
		new_tris[i]->edges[2].twin = &(new_tris[(i + 2) % 3]->edges[1]);
	}
}

DAGNode * insert_pr_on_edge(int pr, int pi_in_tri_index, DAGNode * pi_pj_pk_dagnode, Point * pts) {
	// pr in on edge pi-pj
	Triangle * pi_pj_pk_tri = pi_pj_pk_dagnode->tri;

	int pi = pi_pj_pk_tri->ids[pi_in_tri_index];
	int pj = pi_pj_pk_tri->ids[(pi_in_tri_index + 1) % 3];
	int pk = pi_pj_pk_tri->ids[(pi_in_tri_index + 2) % 3];

	Edge * pj_pi_pl_tri_pj_pi_edge = (&pi_pj_pk_tri->edges[pi_in_tri_index])->twin;
	Triangle * pj_pi_pl_tri = pj_pi_pl_tri_pj_pi_edge->inc;
	DAGNode * pj_pi_pl_dagnode = pj_pi_pl_tri->dag_node;

	int pj_pi_pl_tri_pj_index = pj_pi_pl_tri_pj_pi_edge->start_vertex_index;
	int pl = pj_pi_pl_tri->ids[(pj_pi_pl_tri_pj_index + 2) % 3]; //L

	int p_verteices[4] = { pj, pk, pi, pl };
	DAGNode * dag_belongs_to[4] = { 
		pi_pj_pk_dagnode , 
		pi_pj_pk_dagnode , 
		pj_pi_pl_dagnode , 
		pj_pi_pl_dagnode 
	};
	Edge * outer_twin_edges[4] = { 
		pi_pj_pk_tri->edges[(pi_in_tri_index + 1) % 3].twin,
		pi_pj_pk_tri->edges[(pi_in_tri_index + 2) % 3].twin,
		pj_pi_pl_tri->edges[(pj_pi_pl_tri_pj_index + 1) % 3].twin,
		pj_pi_pl_tri->edges[(pj_pi_pl_tri_pj_index + 2) % 3].twin,
	};
	Triangle * new_tris[4];

	// edge pointers
	for (int i = 0; i < 4; i++) {
		new_tris[i] = new Triangle(p_verteices[i], p_verteices[(i + 1) % 4], pr, pts);
		dag_belongs_to[i]->children[i % 2] = new DAGNode;

		dag_belongs_to[i]->children[i % 2]->tri = new_tris[i];
		new_tris[i]->dag_node = dag_belongs_to[i]->children[i % 2];
		if (outer_twin_edges[i] != nullptr) {
			outer_twin_edges[i]->twin = & new_tris[i]->edges[0];
		}	
		new_tris[i]->edges[0].twin = outer_twin_edges[i];
	}

	for (int i = 0; i < 4; i++) {
		new_tris[i]->edges[1].twin = &(new_tris[(i + 1) % 4]->edges[2]);
		new_tris[i]->edges[2].twin = &(new_tris[(i + 3) % 4]->edges[1]);
	}

	return pj_pi_pl_dagnode;
}

void result_count(DAGNode * head, Point *pts, int n, ll &result, ll &leaf_count) {
	if (head == nullptr) return;
	if (head->visited) return;

	head->visited = true;	

	if (head->is_leaf()) {
		//printf("leaf: (%d,%d,%d)\n", head->tri->ids[0], head->tri->ids[1], head->tri->ids[2]);
		int ids[3];
		for (int i = 0; i < 3; i++) {
			ids[i] = head->tri->ids[i];
		}
		
		for (int i = 0; i < 3; i++) {
			if (ids[i] < n && ids[(i+1) % 3] < n && head->tri->edges[i].counted == false) {
				head->tri->edges[i].counted = true;
				Edge * edge_twin = head->tri->edges[i].twin;
				if (edge_twin != nullptr) { edge_twin->counted = true; }

				leaf_count++;
				result += pts[ids[i]].true_id + pts[ids[(i + 1) % 3]].true_id;
			}
		}
	}
	else {
		for (int i = 0; i < 3; i++) {
			result_count(head->children[i], pts, n, result, leaf_count);
		}
	}
}

void random_increment_delaunay(Point *pts, int n) {
	// find p0
	int y_max = -POINT_MAX;
	int y_max_x_max = -POINT_MAX;
	int y_max_id = -1;
	for (int i = 0; i < n; i++) {
		if (pts[i].y > y_max) {
			y_max_id = i;
			y_max = pts[i].y;
			y_max_x_max = pts[i].x;
		}
		else if (pts[i].y == y_max) {
			if (pts[i].x > y_max_x_max) {
				y_max_id = i;
				y_max_x_max = pts[i].x;
			}
		}
	}

	pts[n].x = POINT_MAX; pts[n].y = -POINT_MAX;
	pts[n+1].x = -POINT_MAX; pts[n+1].y = POINT_MAX;

	swap(pts[y_max_id], pts[n - 1]);
	// random shuffle
	random_shuffle(pts, pts + n - 1);

	//pts[n]:p -1; pts[n]: p -2;
	DAGNode * dag_head = new DAGNode;
	dag_head->tri = new Triangle(n - 1, n + 1, n);
	dag_head->tri->dag_node = dag_head;

	for (int pr = 0; pr < n - 1; pr++) {
		// find triangle contain p_r;
		DAGNode * pi_pj_pk_dagnode = find_triangle_contain_pr(dag_head, pts[pr], pts, n);
		Triangle * pi_pj_pk_tri = pi_pj_pk_dagnode->tri;

		ll to_left_values[3];
		for (int i = 0; i < 3; i++) {
			if (pi_pj_pk_tri->ids[i] < n && pi_pj_pk_tri->ids[(i + 1) % 3] < n) {
				to_left_values[i] = to_left_value(pts[pi_pj_pk_tri->ids[i]], pts[pi_pj_pk_tri->ids[(i + 1) % 3]], pts[pr]);
			}
			else {
				to_left_values[i] = 1;
			}
		}
		if (to_left_values[0] > 0 && to_left_values[1] > 0 && to_left_values[2] > 0) {
			// inside
			insert_pr_inside(pr, pi_pj_pk_dagnode, pts);

			// legalize
			legalize(2, 0, 1, pi_pj_pk_dagnode->children[0], pts, n);
			legalize(2, 0, 1, pi_pj_pk_dagnode->children[1], pts, n);
			legalize(2, 0, 1, pi_pj_pk_dagnode->children[2], pts, n);
		}
		else {
			// on an edge
			int pi_in_tri_index; // edge pi-pj
			for (pi_in_tri_index = 0; pi_in_tri_index < 3; pi_in_tri_index++) {
				if (to_left_values[pi_in_tri_index] == 0) {
					break;
				}
			}
			DAGNode * pj_pi_pl_dagnode = insert_pr_on_edge(pr, pi_in_tri_index, pi_pj_pk_dagnode, pts);
				
			legalize(2, 0, 1, pi_pj_pk_dagnode->children[0], pts, n);
			legalize(2, 0, 1, pi_pj_pk_dagnode->children[1], pts, n);
			legalize(2, 0, 1, pj_pi_pl_dagnode->children[0], pts, n);
			legalize(2, 0, 1, pj_pi_pl_dagnode->children[1], pts, n);
		}
	}


	ll result=0;
	ll leaf_count=0;

	result_count(dag_head, pts, n, result, leaf_count);
	result = result % (leaf_count + 1);
	printf("%lld\n", result);
}

int main() {
	int n;
	scanf("%d", &n);
	Point * pts = new Point[n + 10];
	for (int i = 0; i < n; i++) {
		scanf("%d %d", &pts[i].x, &pts[i].y);
		pts[i].true_id = i + 1;
	}

	random_increment_delaunay(pts, n);
	delete[] pts;
	return 0;
}