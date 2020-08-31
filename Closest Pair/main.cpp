#include <cstdio>
#include <algorithm>
#include <math.h>
using namespace std;

const double MAX_DISTANCE = 1e30;

class Point2D {
public:
	int x;
	int y;
public:
	Point2D() { x = 0; y = 0; }
	Point2D(int x, int y) { x = x; y = y; }
	~Point2D() {}
};

typedef Point2D Point;

bool pts_sort(Point2D &a, Point2D &b) {
	if (a.x < b.x)
		return true;
	else if (a.x > b.x)
		return false;
	else return a.y < b.y;
}

double distance(Point2D a, Point2D b) {
	double tmp = ((double)a.x - b.x) * ((double)a.x - b.x) + ((double)a.y - b.y) * ((double)a.y - b.y);
	return sqrt((double)tmp);
}

double abs_diff(int y1, int y2) {
	double tmp = double(y1) - double(y2);
	return tmp > 0 ? tmp : -tmp;
}

double min_dis(Point2D * pts, int start, int end, Point2D * left, Point2D * right) {
	if (start >= end) {
		return MAX_DISTANCE;
	}

	if (start == end - 1) {
		// important: merge sort according to y
		if (pts[start].y > pts[end].y) {
			swap(pts[start], pts[end]);
		}
		return distance(pts[start], pts[end]);
	}

	int mid = (start + end) / 2;
	int mid_x = pts[mid].x;
	//divide
	double d1 = min_dis(pts, start, mid - 1, left, right);
	double d2 = min_dis(pts, mid, end, left, right);

	double d_min = d1 < d2 ? d1 : d2;

	int left_count = 0;
	int right_count = 0;
	//left&right insert
	for (int i = start; i <= mid - 1; i++) {
		if ((double)mid_x - (double)pts[i].x <= d_min) {
			left[left_count] = pts[i];
			left_count++;
		}
	}
	for (int i = mid; i <= end; i++) {
		if ((double)pts[i].x - (double)mid_x <= d_min) {
			right[right_count] = pts[i];
			right_count++;
		}
	}

	// mid min distance
	int i, j;
	j = 0;
	for (i = 0; i < left_count; i++) {
		int y_i = left[i].y;
		if (j >= right_count)
			break;

		if (abs_diff(y_i, right[j].y) <= d_min) {
			//find end
			int test_start_j = j;
			int test_end_j = test_start_j;
			while (test_end_j < right_count && abs_diff(y_i,right[test_end_j].y) <= d_min) { 
				test_end_j++; 
			}
			for (int k = test_start_j; k < test_end_j; k++) {
				d_min = min(d_min, distance(left[i], right[k]));
			}
		}
		else if(right[j].y < y_i) {
			//find first j
			bool quit_flag = false;
			while (j < right_count) {
				if (abs_diff(y_i, right[j].y) <= d_min) {
					break;
				}
				else if (right[j].y > y_i) { // right[j].y > y_i + d_min
					quit_flag = true;
					break;
				}
				j++;
			}
			if (j >= right_count) {
				break;
			}
			if (quit_flag) {
				continue;
			}

			int test_start_j = j;

			//find end
			int test_end_j = test_start_j;
			while (test_end_j < right_count && abs_diff(y_i,right[test_end_j].y) <= d_min) { 
				test_end_j++; 
			}
			for (int k = test_start_j; k < test_end_j; k++) {
				d_min = min(d_min, distance(left[i], right[k]));
			}
		}
	}

	// merge sort
	// copy
	for (i = start; i <= end; i++) {
		left[i] = pts[i];
	}

	// merge
	i = start;
	j = mid;
	int k = start;
	while (i <= mid-1 && j <= end) {
		if (left[i].y < left[j].y) {
			pts[k] = left[i];
			i++;
			k++;
		}
		else { // left[i].y >= left[j].y
			pts[k] = left[j];
			j++;
			k++;
		}
	}

	for (; i <= mid - 1; i++, k++) {
		pts[k] = left[i];
	}
	for (; j <= end; j++, k++) {
		pts[k] = left[j];
	}

	return d_min;
}

int main() {
	int n = 1;
	scanf("%d", &n);
	Point2D * pts = new Point2D[n];
	Point2D * left = new Point2D[n];
	Point2D * right = new Point2D[n];
	for (int i = 0; i < n; i++) {
		scanf("%d %d", &pts[i].x, &pts[i].y);
	}

	sort(pts, pts + n, pts_sort);
	double d_min = min_dis(pts, 0, n - 1, left, right);
	printf("%.2f", d_min);

	delete[] pts;
	delete[] left;
	delete[] right;
	return 0;
}