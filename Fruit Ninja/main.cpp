#include <cstdio>
#include <algorithm>
#include <vector>
using namespace std;

#define POINT_MAX_X 1e7
#define POINT_MAX_Y 1e14
#define K_INFINITY 1e7
#define RANK_INFINITY 1e30

#define eps 1e-7

#define SAME true
#define OPPOSE false

typedef long long ll;

class Point {
public:
	int x;
	int y;
public:
	Point() { x = 0; y = 0; }
	Point(int mx, int my) { x = mx; y = my; }
	~Point() {}
	const Point& operator =(const Point &b) { x = b.x; y = b.y; return b; }
	bool operator ==(Point &b) { return x == b.x && y == b.y; }
};

class Point_Double {
public:
	ll up_x;
	ll up_y;
	int div;
public:
	Point_Double() { up_x = 0; up_y = 0; div = 1; }
	Point_Double(ll mx, ll my, int mdiv = 1) { up_x = mx; up_y = my; div = mdiv; }
	~Point_Double() {}
	const Point_Double& operator =(const Point_Double &b) {
		up_x = b.up_x; 
		up_y = b.up_y;
		div = b.div;
		return b; 
	}
	bool operator ==(Point_Double &b) { 
		bool x_same = (up_x * (ll)b.div == b.up_x * (ll)div);
		bool y_same = (up_y * (ll)b.div == b.up_y * (ll)div);
		return x_same && y_same;
	}
};

int cross_product(int p1x, int p1y, int p2x, int p2y) {
	return (p1x * p2y) - (p1y * p2x);
}

ll cross_product(ll p1x, ll p1y, ll p2x, ll p2y) {
	return (p1x * p2y) - (p1y * p2x);
}

bool to_left(Point a, Point b, Point c) {
	ll b_a_x = b.x - a.x;
	ll b_a_y = b.y - a.y;
	ll c_a_x = c.x - a.x;
	ll c_a_y = c.y - a.y;
	ll v = b_a_x * c_a_y - b_a_y * c_a_x;
	return v >= 0;
}

bool to_left(Point a, Point b, Point_Double c) {
	ll b_a_x = b.x - a.x;
	ll b_a_y = b.y - a.y;
	ll c_a_x = c.up_x - (ll)a.x * c.div;
	ll c_a_y = c.up_y - (ll)a.y * c.div;
	ll v = b_a_x * c_a_y - b_a_y * c_a_x;
	return v >= 0;
}

class Vertical_Line {
public:
	int x;
	int y1;
	int y2;
};

class Half_Plane {
public:
	Point a;
	Point b;
	int k;
	int m;
public:
	bool contain_point(Point &q) {
		return to_left(a, b, q);
	}
	bool contain_point(Point_Double &q) {
		if (q.up_y != -POINT_MAX_Y)
			return to_left(a, b, q);
		else {
			if (this->a.x < this->b.x) return true;
			else return false;
		}
	}

	Half_Plane() {}
	~Half_Plane() {}

	void assign(Point ma, Point mb,  int mk, int mm) {
		a = ma;
		b = mb;
		k = mk;
		m = mm;
	}
};

bool objective_function_cmp(Point_Double a, Point_Double b) {
	ll ax = a.up_x * b.div;
	ll bx = b.up_x * a.div;
	if (ax > bx) return true;
	else if (ax < bx) return false;
	else {
		ll ay = a.up_y * b.div;
		ll by = b.up_y * a.div;
		return ay < by;
	}
}

double pt_rank(Half_Plane * pa, Point_Double &point) {
	// pa->k != K_INFINITY
	double pointx = (double)point.up_x / (double)point.div;
	return (pa->b.x - pa->a.x) > 0 ? (pointx - pa->a.x) : - (pointx - pa->a.x);
}

void line_intersect(Half_Plane *pa, Half_Plane *pb, bool &valid, bool &has_intersection, Point_Double &intersection, bool &direction) {
	if (pa->k == pb->k) {
		has_intersection = false;
		if (!pa->contain_point(pb->a) && !pb->contain_point(pa->a))
			valid = false;
		else
			valid = true;
		return;
	}
	else {
		valid = true;
		has_intersection = true;
		
		if (pb->k == K_INFINITY) {
			intersection.up_x = pb->a.x;
			intersection.up_y = pa->k * intersection.up_x + pa->m;
			intersection.div = 1;
		}
		else {
			intersection.div = pa->k - pb->k;
			intersection.up_x = pb->m - pa->m;
			intersection.up_y = (ll)pa->k * pb->m - (ll)pb->k * pa->m;

			if (intersection.div < 0) {
				intersection.div = -intersection.div;
				intersection.up_x = -intersection.up_x;
				intersection.up_y = -intersection.up_y;
			}
		}

		int vector_a_x = pa->b.x - pa->a.x;
		int vector_a_y = pa->b.y - pa->a.y;
		int vector_b_x = pb->b.x - pb->a.x;
		int vector_b_y = pb->b.y - pb->a.y;
		int v = cross_product(vector_a_x, vector_a_y, vector_b_x, vector_b_y);
		if (v >= 0) direction = OPPOSE;
		else direction = SAME;
		return;
	}
}

char judge_fruit_ninja(Vertical_Line * lines, int n) {
	Half_Plane * half_planes = new Half_Plane[n * 2 + 10];

	// y1 < y2
	for (int i = 0; i < n; i++) {
		half_planes[i * 2].assign(Point(1, lines[i].x - lines[i].y1), Point(0, -lines[i].y1), lines[i].x , -lines[i].y1);
		half_planes[i * 2 + 1].assign(Point(0, -lines[i].y2), Point(1, lines[i].x - lines[i].y2), lines[i].x, -lines[i].y2);
	}
	half_planes[2 * n].assign(Point(POINT_MAX_X, 0), Point(POINT_MAX_X, 1), K_INFINITY, K_INFINITY);
	//half_planes[2 * n + 1].assign(Point(0, -POINT_MAX_Y), Point(1, -POINT_MAX_Y), 0, -POINT_MAX_Y);

	random_shuffle(half_planes, half_planes + (n * 2));
	random_shuffle(half_planes, half_planes + (n * 2));

	Half_Plane ** included_half_planes = new Half_Plane *[n * 2 + 10];
	Point_Double min_objective_function_point = Point_Double(POINT_MAX_X,-POINT_MAX_Y);
	included_half_planes[0] = &half_planes[2 * n];
	//included_half_planes[1] = &half_planes[2 * n + 1];

	Point_Double intersection;
	bool valid;
	bool has_intersection;
	bool direction;

	const int n_2 = n * 2;
	for (int i = 0; i < n_2; i++) {
		Half_Plane * current_plane = &half_planes[i];
		included_half_planes[i + 1] = current_plane;
		if (current_plane->contain_point(min_objective_function_point))
			continue;
		else {
			double same_largest_rank = -RANK_INFINITY;
			Point_Double same_largest;
			double oppose_smallest_rank = RANK_INFINITY;
			Point_Double oppose_smallest;

			int included_size = i + 1;
			for (int j = 0; j < included_size; j++) {
				line_intersect(current_plane, included_half_planes[j], valid, has_intersection, intersection, direction);

				if (valid == false)
					return 'N';

				// valid == true
				if (!has_intersection) continue;

				double intersection_rank = pt_rank(current_plane, intersection);

				if (direction == SAME) {
					if (same_largest_rank == -RANK_INFINITY || intersection_rank > same_largest_rank) {
						same_largest_rank = intersection_rank;
						same_largest = intersection;
					}
				}
				else {
					// direction == OPPOSE
					if (oppose_smallest_rank == RANK_INFINITY || intersection_rank < oppose_smallest_rank) {
						oppose_smallest_rank = intersection_rank;
						oppose_smallest = intersection;
					}
				}
			}

			if (same_largest_rank > oppose_smallest_rank + eps)
				return 'N';
			else {
				if (same_largest_rank == -RANK_INFINITY)
					min_objective_function_point = oppose_smallest;
				else if (oppose_smallest_rank == RANK_INFINITY)
					min_objective_function_point = same_largest;
				else 
					min_objective_function_point = objective_function_cmp(same_largest, oppose_smallest)
							? same_largest : oppose_smallest;
			}
		}
	}

	delete[] half_planes;
	delete[] included_half_planes;
	return 'Y';
}

int main() {
	int m;
	scanf("%d", &m);
	Vertical_Line * lines; 

	int n;
	char * result = new char[m];
	for (int i = 0; i < m; i++) {
		scanf("%d", &n);
		lines = new Vertical_Line[n];
		for (int j = 0; j < n; j++) {
			scanf("%d %d %d", &lines[j].x, &lines[j].y1, &lines[j].y2);
			if (lines[j].y1 > lines[j].y2) swap(lines[j].y1, lines[j].y2);
		}
		result[i] = judge_fruit_ninja(lines, n);

		delete[] lines;
	}

	for (int i = 0; i < m; i++) {
		putchar(result[i]);
	}
	
	return 0;
}