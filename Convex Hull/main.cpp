#include <cstdio>
#include <algorithm>
using namespace std;

#define MAX_LENGTH 100000
#define XY_MAX 1000000


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
	M_Stack() { data = new Point[100]; capacity = 100; count = 0;  }
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
};

int main() {
	int n;
	scanf("%d", &n);
	Point2D* pts = new Point2D[n+5];
	for (int i = 0; i < n; i++) {
		scanf("%d %d", &pts[i].x, &pts[i].y);
		pts[i].index = i+1;
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

	//reverse first edge
	int first_edge_count = 0;
	int first_edge_x = pts[0].x - p1_x;
	int first_edge_y = pts[0].y - p1_y;
	for (first_edge_count = 1; first_edge_count < n - 1; first_edge_count++) {
		long long cross_p = cross_product(first_edge_x, first_edge_y, pts[first_edge_count].x - p1_x, pts[first_edge_count].y - p1_y);
		if (cross_p != 0)
			break;
	}

	long long M = 1;
	//sort first edge
	if (first_edge_count >= 2) {
		if (first_edge_count == n - 1) {
			M = n;
			for (int i = 1; i <= n; i++) {
				M = (M * (long long)i) % 1000000007;
			}

			printf("%lld\n", M);

			delete[] pts;
			return 0;
		}

		//sort(pts, pts + first_edge_count, first_edge_sort);
		for (int i = 0; i < first_edge_count / 2; i++) {
			swap(pts[i], pts[first_edge_count - i - 1]);
		}

		//for (int i = 1; i < first_edge_count; i++) {
		//	M = (M * (long long)pts[i].index) % 1000000007;
		//}
	}

	// to stack
	M_Stack s(n);
	M_Stack t(n);

	// initialize
	s.push(p1_x, p1_y, p1_index+1);
	s.push(pts[0].x, pts[0].y, pts[0].index);
	for (int i = n - 2 ; i >= 1; i--) {
	//for (int i = n - 2; i >= first_edge_count; i--) {
		t.push(pts[i].x, pts[i].y, pts[i].index);
	}

	while (!t.empty()) {
		Point w = t.top();
		Point v = s.top();
		Point u = s.second_top();
		if (to_left(u, v, w)) {
			s.push(t.pop());
		}
		else {
			s.pop();
		}
	}

	int h = s.size();
	//int h = s.size() + first_edge_count - 1;
	M = (M * (long long)h) % 1000000007;
	while (!s.empty()) {
		M = (M * (long long)(s.top().index)) % 1000000007;
		s.pop();
	}

	printf("%lld\n", M);

	delete [] pts;
	return 0;
}