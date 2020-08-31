#include <cstdio>
#include <algorithm>
#include <vector>
using namespace std;

#define TRAPEZOID_NODE 0
#define X_NODE 1
#define Y_NODE 2
#define POINT_MAX 2000000

typedef long long ll;


class Point2D {
public:
	int x;
	int y;
	int true_id;
public:
	Point2D() { x = 0; y = 0; }
	Point2D(int mx, int my) { x = mx; y = my; }
	~Point2D() {}
	const Point2D& operator =(const Point2D &b) { x = b.x; y = b.y; return b; }
	bool operator ==(Point2D &b) { return x == b.x && y == b.y; }
};

typedef Point2D Point;

class Edge {
public:
	Point *left;
	Point *right;
	int true_id;
public:
	Edge() { left = nullptr; right = nullptr; true_id = -1; }
	Edge(Point *l, Point *r) { left = l; right = r; true_id = -1; }
	~Edge() {}
	const Edge& operator =(const Edge &b) { left = b.left; right = b.right; true_id = b.true_id; return b; }
};

bool point_smaller(Point a, Point b) {
	if (a.x < b.x) return true;
	else {
		if (a.x > b.x) return false;
		else return a.y < b.y;
	}
}

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

bool to_left(Point a, Edge edge) {
	ll v = to_left_value(*edge.left, *edge.right, a);
	return v > 0;
}

class Map_Node;
typedef Map_Node Trapezoid_Node;
typedef Map_Node X_Node;
typedef Map_Node Y_Node;

class Map_Node {
public:
	int type;
	Map_Node * left_child = nullptr;
	Map_Node * right_child = nullptr;
public:
	// X_Node
	Point * pt = nullptr;
	void assign_X_Node(Point *p) { type = X_NODE; pt = p; }
public:
	// Y_Node
	Edge * edge = nullptr;
	void assign_Y_Node(Edge *e) { type = Y_NODE; edge = e; }
public:
	// Trapezoid_Node
	Edge * top = nullptr;
	Edge * bottom = nullptr;
	Point * leftp = nullptr;
	Point * rightp = nullptr;

	Trapezoid_Node * left_neighbors[2] = { nullptr,nullptr };
	Trapezoid_Node * right_neighbors[2] = { nullptr,nullptr };
	void assign_Trapezoid_Node() { type = TRAPEZOID_NODE; }
	void assign_Trapezoid_Node(Edge * n_top, Edge * n_bottom, Point * n_leftp, Point *n_rightp) {
		type = TRAPEZOID_NODE;
		top = n_top;
		bottom = n_bottom;
		leftp = n_leftp;
		rightp = n_rightp;
	}
	void assign_Trapezoid_Node(Trapezoid_Node &m) {
		type = TRAPEZOID_NODE;
		top = m.top;
		bottom = m.bottom;
		leftp = m.leftp;
		rightp = m.rightp;
		copy_left_neighbor(m);
		copy_right_neighbor(m);
	}
	void copy_left_neighbor(Trapezoid_Node &m) {
		left_neighbors[0] = m.left_neighbors[0];
		left_neighbors[1] = m.left_neighbors[1];
	}
	void copy_right_neighbor(Trapezoid_Node &m) {
		right_neighbors[0] = m.right_neighbors[0];
		right_neighbors[1] = m.right_neighbors[1];
	}
public:
	Map_Node() { type = -1; left_child = nullptr; right_child = nullptr; }
	~Map_Node() {}
	bool is_leaf() { return type == TRAPEZOID_NODE; }
	void link_left_child(Map_Node * a) { left_child = a; }
	void link_right_child(Map_Node * a) { right_child = a; }
	bool smaller(Point query) { 
		if (type == X_NODE)
			return point_smaller(query, *pt);
		else
			// type == Y_NODE
			return to_left(query, *edge);
	}
};

Map_Node * strict_query_point(Point p, Map_Node * head) {
	Map_Node * node_pt = head;
	while (!node_pt->is_leaf()) {
		if (node_pt->type == X_NODE) {
			// pt->type = X_NODE
			if (*node_pt->pt == p) // at a point
				break;

			if (node_pt->smaller(p)) node_pt = node_pt->left_child;
			else node_pt = node_pt->right_child;
		}
		else {
			// pt->type = Y_NODE
			Y_Node * m = (Y_Node *)node_pt;
			ll v = to_left_value(*m->edge->left, *m->edge->right, p);

			if (v > 0)  node_pt = node_pt->left_child;
			else if (v < 0) node_pt = node_pt->right_child;
			else // on a edge
				break;
		}
	}

	return node_pt;
}

Trapezoid_Node * query_point(Point * p, Map_Node * head) {
	Map_Node * pt = head;
	while (!pt->is_leaf()) {
		if (pt->smaller(*p)) pt = pt->left_child;
		else pt = pt->right_child;
	}

	return (Trapezoid_Node *)pt;
}

void find_all_trapezoids(Edge * new_edge, Map_Node * head, vector<Trapezoid_Node *> &trapezoids, Trapezoid_Node * &left_point_trapezoid, Trapezoid_Node * &right_point_trapezoid) {
	trapezoids.clear();

	Point * new_edge_left_point = new_edge->left;   //p
	Point * new_edge_right_point = new_edge->right; //q
	Map_Node * left_point_trapezoid_parent;
	left_point_trapezoid = query_point(new_edge_left_point, head);
	trapezoids.push_back(left_point_trapezoid);

	Trapezoid_Node * trapezoid_pointer = left_point_trapezoid;
	while (point_smaller(*trapezoid_pointer->rightp, *new_edge_right_point)) {
		if (to_left(*trapezoid_pointer->rightp, *new_edge))
			trapezoid_pointer = trapezoid_pointer->right_neighbors[1];
		else
			trapezoid_pointer = trapezoid_pointer->right_neighbors[0];

		if (trapezoid_pointer == nullptr)
			break;
		else
			trapezoids.push_back(trapezoid_pointer);
	}
	
	right_point_trapezoid = trapezoids.back();
}

void set_single_left_neighbor_neighbor_pointer(int index, Trapezoid_Node * current_trapezoid, Trapezoid_Node * new_trapezoid, Trapezoid_Node * banned = nullptr) {
	new_trapezoid->left_neighbors[index] = current_trapezoid->left_neighbors[index];
	if (current_trapezoid->left_neighbors[index] != nullptr && current_trapezoid->left_neighbors[index] != banned) {
		Trapezoid_Node * p = current_trapezoid->left_neighbors[index];
		for (int j = 0; j < 2; j++) {
			if (p->right_neighbors[j] == current_trapezoid)
				p->right_neighbors[j] = new_trapezoid;
		}
	}
}

void set_single_right_neighbor_neighbor_pointer(int index, Trapezoid_Node * current_trapezoid, Trapezoid_Node * new_trapezoid, Trapezoid_Node * banned = nullptr) {
	new_trapezoid->right_neighbors[index] = current_trapezoid->right_neighbors[index];
	if (current_trapezoid->right_neighbors[index] != nullptr && current_trapezoid->right_neighbors[index] != banned) {
		Trapezoid_Node * p = current_trapezoid->right_neighbors[index];
		for (int j = 0; j < 2; j++) {
			if (p->left_neighbors[j] == current_trapezoid)
				p->left_neighbors[j] = new_trapezoid;
		}
	}
}

void set_left_neighbors_neighbor_pointer(Trapezoid_Node * current_trapezoid, Trapezoid_Node * new_trapezoid, Trapezoid_Node * banned=nullptr) {
	set_single_left_neighbor_neighbor_pointer(0, current_trapezoid, new_trapezoid, banned);
	set_single_left_neighbor_neighbor_pointer(1, current_trapezoid, new_trapezoid, banned);
}

void set_right_neighbors_neighbor_pointer(Trapezoid_Node * current_trapezoid, Trapezoid_Node * new_trapezoid, Trapezoid_Node * banned = nullptr) {
	set_single_right_neighbor_neighbor_pointer(0, current_trapezoid, new_trapezoid, banned);
	set_single_right_neighbor_neighbor_pointer(1, current_trapezoid, new_trapezoid, banned);
}

void set_toporbottom_neighbors_neighbor_pointer(int index, Trapezoid_Node * current_trapezoid, Trapezoid_Node * new_trapezoid, Trapezoid_Node * banned_left = nullptr, Trapezoid_Node * banned_right = nullptr) {
	set_single_left_neighbor_neighbor_pointer(index, current_trapezoid, new_trapezoid, banned_left);
	set_single_right_neighbor_neighbor_pointer(index, current_trapezoid, new_trapezoid, banned_right);
}

Map_Node * build_trapezoid_map(Point * pts, Edge * edges, int n) {
	Map_Node * head = nullptr;

	random_shuffle(edges, edges + n);

	// create R
	pts[2 * n] = Point(-POINT_MAX, POINT_MAX);
	pts[2 * n + 1] = Point(POINT_MAX, POINT_MAX);
	pts[2 * n + 2] = Point(-POINT_MAX, -POINT_MAX);
	pts[2 * n + 3] = Point(-POINT_MAX, POINT_MAX);
	edges[n] = Edge(&pts[2 * n], &pts[2 * n + 1]);
	edges[n + 1] = Edge(&pts[2 * n + 2], &pts[2 * n + 3]);

	head = new Trapezoid_Node();
	head->assign_Trapezoid_Node(&edges[n], &edges[n + 1], &pts[2 * n], &pts[2 * n + 1]);

	Trapezoid_Node * left_point_trapezoid;
	Trapezoid_Node * right_point_trapezoid;
	vector<Trapezoid_Node *> trapezoids;
	for (int i = 0; i < n; i++) {
		find_all_trapezoids(&edges[i], head, trapezoids, left_point_trapezoid, right_point_trapezoid);
		if (trapezoids.size() == 1) {
			// 2 end points in 1 trapezoid
			{
				Edge * current_edge = &edges[i];
				Trapezoid_Node * current_trapezoid = trapezoids[0];
				Trapezoid_Node * A = new Trapezoid_Node; A->assign_Trapezoid_Node(
					current_trapezoid->top, current_trapezoid->bottom,
					current_trapezoid->leftp, current_edge->left
				);
				set_left_neighbors_neighbor_pointer(current_trapezoid, A);

				Trapezoid_Node * B = new Trapezoid_Node; B->assign_Trapezoid_Node(
					current_trapezoid->top, current_trapezoid->bottom,
					current_edge->right, current_trapezoid->rightp
				);
				set_right_neighbors_neighbor_pointer(current_trapezoid, B);

				Trapezoid_Node * C = new Trapezoid_Node; C->assign_Trapezoid_Node(
					current_trapezoid->top, current_edge,
					current_edge->left, current_edge->right
				);

				Trapezoid_Node * D = new Trapezoid_Node; D->assign_Trapezoid_Node(
					current_edge, current_trapezoid->bottom,
					current_edge->left, current_edge->right
				);

				A->right_neighbors[0] = C;
				A->right_neighbors[1] = D;
				B->left_neighbors[0] = C;
				B->left_neighbors[1] = D;
				C->left_neighbors[0] = A;
				C->right_neighbors[0] = B;
				D->left_neighbors[1] = A;
				D->right_neighbors[1] = B;

				// root == p
				X_Node * p = current_trapezoid; p->assign_X_Node(current_edge->left);
				X_Node * q = new X_Node;		q->assign_X_Node(current_edge->right);
				Y_Node * si = new Y_Node;		si->assign_Y_Node(current_edge);

				p->link_left_child(A);
				p->link_right_child(q);
				q->link_left_child(si);
				q->link_right_child(B);
				si->link_left_child(C);
				si->link_right_child(D);
			}
		}
		else {
			// multiple trapezoids
			{
				int trapezoids_size = trapezoids.size();
				Edge * current_edge = &edges[i];
				// deal with left most
				Trapezoid_Node * A = new Trapezoid_Node; A->assign_Trapezoid_Node(
					left_point_trapezoid->top, left_point_trapezoid->bottom,
					left_point_trapezoid->leftp, current_edge->left
				);
				set_left_neighbors_neighbor_pointer(left_point_trapezoid, A);

				Trapezoid_Node * top_start = new Trapezoid_Node; top_start->assign_Trapezoid_Node(
					left_point_trapezoid->top, current_edge,
					current_edge->left, left_point_trapezoid->rightp
				);
				set_single_right_neighbor_neighbor_pointer(0, left_point_trapezoid, top_start, trapezoids[1]);

				Trapezoid_Node * bottom_start = new Trapezoid_Node; bottom_start->assign_Trapezoid_Node(
					current_edge, left_point_trapezoid->bottom,
					current_edge->left, left_point_trapezoid->rightp
				);
				set_single_right_neighbor_neighbor_pointer(1, left_point_trapezoid, bottom_start, trapezoids[1]);

				vector<Trapezoid_Node *> top_trapezoids;
				vector<Trapezoid_Node *> distinct_top_trapezoids;
				vector<Trapezoid_Node *> bottom_trapezoids;
				vector<Trapezoid_Node *> distinct_bottom_trapezoids;
				top_trapezoids.push_back(top_start);
				bottom_trapezoids.push_back(bottom_start);
				distinct_top_trapezoids.push_back(top_start);
				distinct_bottom_trapezoids.push_back(bottom_start);

				for (int j = 1; j < trapezoids_size; j++) {
					Trapezoid_Node * current_trapezoid = trapezoids[j];
					Trapezoid_Node * current_top = top_trapezoids.back();
					Trapezoid_Node * current_bottom = bottom_trapezoids.back();

					// top
					if (current_trapezoid->top == current_top->top) {
						// top merge
						if (j < trapezoids_size - 1) {
							current_top->rightp = current_trapezoid->rightp;
							set_single_right_neighbor_neighbor_pointer(0, current_trapezoid, current_top, trapezoids[j + 1]);
						}
						else {
							// j == trapezoid_size - 1
							current_top->rightp = current_edge->right;
						}
						top_trapezoids.push_back(current_top);
					}
					else {
						// new top
						Trapezoid_Node * new_top = new Trapezoid_Node;
						if (j < trapezoids_size - 1) {
							new_top->assign_Trapezoid_Node(
								current_trapezoid->top, current_edge,
								current_trapezoid->leftp, current_trapezoid->rightp
							);
							set_single_left_neighbor_neighbor_pointer(0, current_trapezoid, new_top, trapezoids[j - 1]);
							set_single_right_neighbor_neighbor_pointer(0, current_trapezoid, new_top, trapezoids[j + 1]);
						}
						else {
							// j == trapezoids_size - 1
							new_top->assign_Trapezoid_Node(
								current_trapezoid->top, current_edge,
								current_trapezoid->leftp, current_edge->right
							);
							set_single_left_neighbor_neighbor_pointer(0, current_trapezoid, new_top, trapezoids[j - 1]);
						}
						top_trapezoids.push_back(new_top);
						distinct_top_trapezoids.push_back(new_top);
					}

					// bottom
					if (current_trapezoid->bottom == current_bottom->bottom) {
						// bottom merge
						if (j < trapezoids_size - 1) {
							current_bottom->rightp = current_trapezoid->rightp;
							set_single_right_neighbor_neighbor_pointer(1, current_trapezoid, current_bottom, trapezoids[j + 1]);
						}
						else {
							// j == trapezoid_size - 1
							current_bottom->rightp = current_edge->right;
						}
						bottom_trapezoids.push_back(current_bottom);
					}
					else {
						// new bottom
						Trapezoid_Node * new_bottom = new Trapezoid_Node; 
						if (j < trapezoids_size - 1) {
							new_bottom->assign_Trapezoid_Node(
								current_edge, current_trapezoid->bottom,
								current_trapezoid->leftp, current_trapezoid->rightp
							);
							set_single_left_neighbor_neighbor_pointer(1, current_trapezoid, new_bottom, trapezoids[j - 1]);
							set_single_right_neighbor_neighbor_pointer(1, current_trapezoid, new_bottom, trapezoids[j + 1]);
						}
						else {
							// j == trapezoid_size - 1
							new_bottom->assign_Trapezoid_Node(
								current_edge, current_trapezoid->bottom,
								current_trapezoid->leftp, current_edge->right
							);
							set_single_left_neighbor_neighbor_pointer(1, current_trapezoid, new_bottom, trapezoids[j - 1]);
						}
						bottom_trapezoids.push_back(new_bottom);
						distinct_bottom_trapezoids.push_back(new_bottom);
					}
				}

				// deal with right most
				Trapezoid_Node * B = new Trapezoid_Node; B->assign_Trapezoid_Node(
					right_point_trapezoid->top, right_point_trapezoid->bottom,
					current_edge->right, right_point_trapezoid->rightp
				);
				set_right_neighbors_neighbor_pointer(right_point_trapezoid, B);

				// compute neighbor pointer for new trapezoids
				A->right_neighbors[0] = distinct_top_trapezoids[0];
				distinct_top_trapezoids[0]->left_neighbors[0] = A;
				A->right_neighbors[1] = distinct_bottom_trapezoids[0];
				distinct_bottom_trapezoids[0]->left_neighbors[1] = A;

				for (int j = 0; j < distinct_top_trapezoids.size() - 1; j++) {
					distinct_top_trapezoids[j]->right_neighbors[1] = distinct_top_trapezoids[j + 1];
					distinct_top_trapezoids[j + 1]->left_neighbors[1] = distinct_top_trapezoids[j];
				}

				for (int j = 0; j < distinct_bottom_trapezoids.size() - 1; j++) {
					distinct_bottom_trapezoids[j]->right_neighbors[0] = distinct_bottom_trapezoids[j + 1];
					distinct_bottom_trapezoids[j + 1]->left_neighbors[0] = distinct_bottom_trapezoids[j];
				}

				distinct_top_trapezoids.back()->right_neighbors[0] = B;
				B->left_neighbors[0] = distinct_top_trapezoids.back();
				distinct_bottom_trapezoids.back()->right_neighbors[1] = B;
				B->left_neighbors[1] = distinct_bottom_trapezoids.back();
			
				// delete space & link new Map_Node s
				for (int j = 0; j < trapezoids_size; j++) {
					if (j == 0) {
						// first trapezoid
						X_Node * p = trapezoids[j]; p->assign_X_Node(current_edge->left);
						Y_Node * si = new Y_Node; si->assign_Y_Node(current_edge);
						
						p->link_left_child(A);
						p->link_right_child(si);
						si->link_left_child(top_trapezoids[j]);
						si->link_right_child(bottom_trapezoids[j]);
					}
					else if (j == trapezoids_size - 1) {
						X_Node * q = trapezoids[j]; q->assign_X_Node(current_edge->right);
						Y_Node * si = new Y_Node; si->assign_Y_Node(current_edge);

						q->left_child = si;
						q->right_child = B;
						si->left_child = top_trapezoids[j];
						si->right_child = bottom_trapezoids[j];
					}
					else {
						Y_Node * si = trapezoids[j]; si->assign_Y_Node(current_edge);
						si->left_child = top_trapezoids[j];
						si->right_child = bottom_trapezoids[j];
					}
				}
			}
		}
		
	}

	return head;
}

int main() {

	int n;
	int m;

	scanf("%d %d", &n, &m);
	Point * pts = new Point[2 * n + 20];
	Edge * edges = new Edge[n + 10];
	for (int i = 0; i < n; i++) {
		scanf("%d %d %d %d", &pts[i * 2].x, &pts[i * 2].y, &pts[i * 2 + 1].x, &pts[i * 2 + 1].y);
		pts[i * 2].true_id = i + 1;
		pts[i * 2 + 1].true_id = i + 1;
		if (point_smaller(pts[i * 2], pts[i * 2 + 1])) {
			edges[i].left = &pts[i * 2];
			edges[i].right = &pts[i * 2 + 1];
		}
		else {
			edges[i].left = &pts[i * 2 + 1];
			edges[i].right = &pts[i * 2];
		}
		edges[i].true_id = i + 1;
	}

	Map_Node * head = build_trapezoid_map(pts, edges, n);

	Point query_point;
	Map_Node *p;
	int * result = new int[m];
	for (int i = 0; i < m; i++) {
		scanf("%d %d", &query_point.x, &query_point.y);
		p = strict_query_point(query_point, head);
		if (p->type == X_NODE) {
			X_Node * tmp = (X_Node *)p;
			result[i] = tmp->pt->true_id;
		}
		else if (p->type == Y_NODE) {
			Y_Node * tmp = (Y_Node *)p;
			result[i] = tmp->edge->true_id;
		}
		else{
			// p->type == TRAPEZOID_NODE
			Trapezoid_Node * tmp = (Trapezoid_Node *)p;
			if (tmp->rightp->x == query_point.x)
				result[i] = tmp->rightp->true_id;
			else 
				result[i] = tmp->top->true_id;
		}
	}

	for (int i = 0; i < m; i++) {
		if (result[i] != -1)printf("%d\n", result[i]);
		else printf("N\n");
	}
	return 0;
}