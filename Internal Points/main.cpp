#include <cstdio>
#include <algorithm>
#include <math.h>
using namespace std;

#define XY_MAX 2000000000

class Point2D {
public:
	int x;
	int y;
	int index;
};

typedef Point2D Point;

int p1_x, p1_y;

long long cross_product(int p1x, int p1y, int p2x, int p2y) {
	return ((long long)p1x * p2y) - ((long long)p1y * p2x);
}

bool cross_sort(Point2D &a, Point2D &b) {
	long long tmp = cross_product(a.x - p1_x, a.y - p1_y, b.x - p1_x, b.y - p1_y);
	if (tmp > 0) {
		return true;
	}
	else if (tmp < 0) {
		return false;
	}
	else {
		long long dis1 = ((long long)a.x - p1_x) * ((long long)a.x - p1_x) + ((long long)a.y - p1_y) * ((long long)a.y - p1_y);
		long long dis2 = ((long long)b.x - p1_x) * ((long long)b.x - p1_x) + ((long long)b.y - p1_y) * ((long long)b.y - p1_y);
		return dis1 > dis2;
	}
}

bool first_edge_sort(Point2D &a, Point2D &b) {
	long long dis1 = ((long long)a.x - p1_x) * ((long long)a.x - p1_x) + ((long long)a.y - p1_y) * ((long long)a.y - p1_y);
	long long dis2 = ((long long)b.x - p1_x) * ((long long)b.x - p1_x) + ((long long)b.y - p1_y) * ((long long)b.y - p1_y);
	return dis1 < dis2;
}

bool to_left(Point &p, Point &q, Point &s) {
	//long long area = (long long)p.x * q.y - (long long)p.y * q.x 
	//	+ (long long)q.x * s.y - (long long)q.y * s.x 
	//	+ (long long)s.x * p.y - (long long)s.y * p.x;

	long long area = cross_product(q.x - p.x, q.y - p.y, s.x - p.x, s.y - p.y);
	return area >= 0;
}


class M_Stack {
public:
	Point * data;
	int count;
	int capacity;
public:
	M_Stack() { data = new Point[100]; capacity = 100; count = 0; }
	M_Stack(int n) { data = new Point[n]; capacity = n; count = 0; }
	~M_Stack() { delete[] data; return; }
	void push(int x, int y, int index) {
		this->data[this->count].x = x;
		this->data[this->count].y = y;
		this->data[this->count].index = index;
		this->count++;
	}
	void push(Point a) {
		this->data[this->count] = a;
		this->count++;
	}
	Point pop() {
		if (this->count <= 0)
			return Point();
		this->count--;
		return this->data[count];
	}
	int size() {
		return this->count;
	}
	Point top() {
		if (this->count <= 0) {
			return Point();
		}
		return this->data[count - 1];
	}
	Point second_top() {
		if (this->count <= 1) {
			return Point();
		}
		return this->data[count - 2];
	}
	bool empty() {
		return this->count == 0;
	}

	//index related
	int next(int i) {
		i = i + 1;
		if (i >= this->count) {
			i = i - this->count;
		}
		return i;
	}

	int pre(int i) {
		i = i - 1;
		if (i < 0) {
			i = i + this->count;
		}
		return i;
	}
};


void graham_scan(M_Stack &s, Point2D * pts, int n, Point2D * inner_pts, int &inner_n , bool record_inner = true) {
	if (n < 3) {
		for (int i = 0; i < n; i++) {
			s.push(pts[i]);
		}
		return;
	}
	
	// find y_min which has the smallest x
	int y_min = XY_MAX;
	int x_min = XY_MAX;
	int p1_index = -1;

	for (int i = 0; i < n; i++) {
		if (pts[i].y < y_min) {
			y_min = pts[i].y;
			x_min = pts[i].x;
			p1_index = i;
		}
		else if (pts[i].y == y_min) {
			if (pts[i].x < x_min) {
				x_min = pts[i].x;
				p1_index = i;
			}
		}
	}

	// pre sorting
	p1_x = pts[p1_index].x;
	p1_y = pts[p1_index].y;

	swap(pts[p1_index], pts[n - 1]);
	sort(pts, pts + n - 1, cross_sort);


	// reverse first edge
	int first_edge_count = 0;
	int first_edge_x = pts[0].x - p1_x;
	int first_edge_y = pts[0].y - p1_y;
	for (first_edge_count = 1; first_edge_count < n - 1; first_edge_count++) {
		long long cross_p = cross_product(first_edge_x, first_edge_y, pts[first_edge_count].x - p1_x, pts[first_edge_count].y - p1_y);
		if (cross_p != 0)
			break;
	}

	//sort first edge
	if (first_edge_count >= 2) {
		if (first_edge_count == n - 1) {
			for (int i = 0; i < n; i++) {
				s.push(pts[i]);
			}
			return;
		}

		for (int i = 0; i < first_edge_count / 2; i++) {
			swap(pts[i], pts[first_edge_count - i - 1]);
		}
	}


	// to stack
	//M_Stack s(n);
	M_Stack t(n);

	// initialize
	s.push(p1_x, p1_y, p1_index + 1);
	s.push(pts[0].x, pts[0].y, pts[0].index);
	for (int i = n - 2; i >= 1; i--) {
		t.push(pts[i].x, pts[i].y, pts[i].index);
	}

	int inner_pts_count = 0;

	while (!t.empty()) {
		Point w = t.top();
		Point v = s.top();
		Point u = s.second_top();
		if (to_left(u, v, w)) {
			s.push(t.pop());
		}
		else {
			// inner pts insert into
			if (record_inner) {
				inner_pts[inner_pts_count] = s.top();
				inner_pts_count++;
			}
			s.pop();
		}
	}

	if (record_inner) {
		inner_n = inner_pts_count;
	}
	return;
}

double dis(Point &a, Point &b, Point &p) {
	int ex = b.x - a.x;
	int ey = b.y - a.y;
	int px = p.x - a.x;
	int py = p.y - a.y;
	
	long long c1 = cross_product(ex, ey, px, py);
	long long e_mol_sq = (long long)ex * ex + (long long)ey * ey;
	double e_mol = sqrt(e_mol_sq);
	return (double)c1 / e_mol;
}

int main() {
	int n;
	scanf("%d", &n);
	Point2D* pts = new Point2D[n + 1];
	Point2D* inner_pts = new Point2D[n + 1];
	for (int i = 0; i < n; i++) {
		scanf("%d %d", &pts[i].x, &pts[i].y);
		pts[i].index = i + 1;
	}

	int inner_n;
	int tmp;
	M_Stack outer_hull(n);
	graham_scan(outer_hull, pts, n, inner_pts, inner_n, true);
	M_Stack inner_hull(inner_n);
	graham_scan(inner_hull, inner_pts, inner_n, NULL, tmp, false);

	printf("%d\n", outer_hull.size());

	
	int outer_hull_count = outer_hull.size();
	int inner_hull_count = inner_hull.size();

	double * results = new double[outer_hull_count];

	if (inner_hull_count <= 4) {
		for (int i = 0; i < outer_hull_count; i++) {
			double min_dis = XY_MAX;
			int i_next = outer_hull.next(i);
			for (int j = 0; j < inner_hull_count; j++) {
				double c_dis = dis(outer_hull.data[i], outer_hull.data[i_next], inner_hull.data[j]);
				if (c_dis < min_dis) {
					min_dis = c_dis;
				}
			}
			results[i] = min_dis;
		}
	}
	else {
		int outer_pointer = 0;
		int inner_pointer = 0;

		int outer_pointer_next;
		int inner_pointer_pre;
		int inner_pointer_next;
		for (outer_pointer = 0; outer_pointer < outer_hull_count; outer_pointer++) {
			outer_pointer_next = outer_hull.next(outer_pointer);

			int ex = outer_hull.data[outer_pointer_next].x - outer_hull.data[outer_pointer].x;
			int ey = outer_hull.data[outer_pointer_next].y - outer_hull.data[outer_pointer].y;

			while (1) {
				inner_pointer_pre = inner_hull.pre(inner_pointer);
				inner_pointer_next = inner_hull.next(inner_pointer);

				int f1x = inner_hull.data[inner_pointer_next].x - inner_hull.data[inner_pointer].x;
				int f1y = inner_hull.data[inner_pointer_next].y - inner_hull.data[inner_pointer].y;

				int f2x = inner_hull.data[inner_pointer_pre].x - inner_hull.data[inner_pointer].x;
				int f2y = inner_hull.data[inner_pointer_pre].y - inner_hull.data[inner_pointer].y;

				long long c1 = cross_product(ex, ey, f1x, f1y);
				long long c2 = cross_product(ex, ey, f2x, f2y);
				if (c1 >= 0 && c2 >= 0) {
					break;
				}
				else if (c1 < 0) {
					inner_pointer = inner_pointer_next;
				}
				else { // c2 < 0
					inner_pointer = inner_pointer_pre;
				}
			}
			results[outer_pointer] = dis(outer_hull.data[outer_pointer], outer_hull.data[outer_pointer_next], inner_hull.data[inner_pointer]);
		}
	}


	sort(results, results + outer_hull_count);
	for (int i = 0; i < outer_hull_count; i++) {
		printf("%lf\n", results[i]);
	}

	delete[] pts;
	delete[] inner_pts;
	return 0;
}