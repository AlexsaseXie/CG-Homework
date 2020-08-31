#include <cstdio>
#include <algorithm>
#include <vector>
using namespace std;

// can use INT_MAX and INT_MIN

#define INT_MAX 2147483647
#define INT_MIN (-INT_MAX-1)

#define NEG_INFINITY false
#define POS_INFINITY true

#define LEFT true
#define RIGHT false

class Point {
public:
	int x;
	int y;
public:
	Point() :x(0), y(0) {}
	Point(int a, int b) { x = a; y = b; }
	~Point() {}
	const Point& operator =(const Point& b) {
		x = b.x;
		y = b.y;
		return b;
	}
};

class Point_And_Link {
public:
	Point * pts;
	int left_child_index;
	int right_child_index;
public:
	Point_And_Link() { pts = nullptr; left_child_index = -1; right_child_index = -1; }
	Point_And_Link(Point * mpts, int mleft_child_index = -1, int mright_child_index = -1) {
		pts = mpts;
		left_child_index = mleft_child_index;
		right_child_index = mright_child_index;
	}
};

class RangeTree_Node {
public:
	Point * split;
	RangeTree_Node *left_child = nullptr;
	RangeTree_Node *right_child = nullptr;
	vector<Point_And_Link> points_and_link;		// sorted by y
public:
	RangeTree_Node() { left_child = nullptr; right_child = nullptr; }
	RangeTree_Node(Point * mx) { split = mx; left_child = nullptr; right_child = nullptr; }
	bool is_leaf() { return left_child == nullptr && right_child == nullptr; }
};

class Window {
public:
	int xmin;
	int xmax;
	int ymin;
	int ymax;
public:
	Window() {}
	~Window() {}
	Window(int m_xmin, int m_xmax, int m_ymin, int m_ymax) { 
		xmin = m_xmin; 
		xmax = m_xmax; 
		ymin = m_ymin; 
		ymax = m_ymax; 
	}
};

bool pts_smaller(Point &a, Point &b) {
	if (a.x < b.x) return true;
	else if (a.x > b.x) return false;
	else return a.y < b.y;
}

class RangeTree {
public:
	RangeTree_Node * head = nullptr;
public:
	RangeTree() { head = nullptr; }
	~RangeTree() {}

	// compare functions
	bool pts_x_smaller(Point &a, int x, bool y_flag) {
		if (a.x < x) return true;
		else if (a.x > x) return false;
		else {
			if (y_flag == NEG_INFINITY) return false;
			else return true;
		}
	}
	bool pts_x_smaller(Point &a, Point &b) {
		if (a.x < b.x) return true;
		else if (a.x > b.x) return false;
		else return a.y < b.y;
	}
	bool pts_y_smaller(Point &a, int y, bool x_flag) {
		if (a.y < y) return true;
		else if (a.y > y) return false;
		else {
			if (x_flag == NEG_INFINITY) return false;
			else return true;
		}
	}
	bool pts_y_smaller(Point &a, Point &b) {
		if (a.y < b.y) return true;
		else if (a.y > b.y) return false;
		else return a.x < b.x;
	}

	// build range tree functions
	RangeTree_Node * build_range_tree(Point * pts, int left, int right) {
		if (left == right) {
			RangeTree_Node * p = new RangeTree_Node;
			p->split = &pts[left];
			p->left_child = nullptr;
			p->right_child = nullptr;
			p->points_and_link.push_back(Point_And_Link(&pts[left]));
			return p;
		}

		int mid = (right + left) / 2;
		RangeTree_Node * p_left_child = build_range_tree(pts, left, mid);
		RangeTree_Node * p_right_child = build_range_tree(pts, mid + 1, right);

		RangeTree_Node * p = new RangeTree_Node;
		p->split = &pts[mid];
		p->left_child = p_left_child;
		p->right_child = p_right_child;

		int left_point_size = p_left_child->points_and_link.size();
		int right_point_size = p_right_child->points_and_link.size();

		int l_pointer = 0;
		int r_pointer = 0;

		while (1) {
			if (l_pointer >= left_point_size || r_pointer >= right_point_size)
				break;

			if (pts_y_smaller(*p_left_child->points_and_link[l_pointer].pts, *p_right_child->points_and_link[r_pointer].pts)) {
				p->points_and_link.push_back(Point_And_Link(p_left_child->points_and_link[l_pointer].pts, l_pointer, r_pointer));
				l_pointer++;
			}
			else {
				p->points_and_link.push_back(Point_And_Link(p_right_child->points_and_link[r_pointer].pts, l_pointer, r_pointer));
				r_pointer++;
			}
		}

		while (l_pointer < left_point_size) {
			p->points_and_link.push_back(Point_And_Link(p_left_child->points_and_link[l_pointer].pts, l_pointer, r_pointer));
			l_pointer++;
		}

		while (r_pointer < right_point_size) {
			p->points_and_link.push_back(Point_And_Link(p_right_child->points_and_link[r_pointer].pts, l_pointer, r_pointer));
			r_pointer++;
		}

		return p;
	}
	void build_range_tree_from_pts(Point * pts, int n, bool need_sort=true) {
		if (need_sort)
			sort(pts, pts + n, pts_smaller);
		this->head = build_range_tree(pts, 0, n - 1);
	}
public:
	// search functions
	int binary_search_y_first_gt(vector<Point_And_Link> &points_and_link, int y, bool x_flag, int n) {
		// return [0, n-1] 
		// or n if query is larger than any point in the vector 

		int low = 0; int high = n - 1; int middle = 0;
		while (low < high)
		{
			middle = (low + high) / 2;
			if (!pts_y_smaller(*points_and_link[middle].pts, y, x_flag)) high = middle;
			else low = middle + 1;
		}

		if (pts_y_smaller(*points_and_link[high].pts, y, x_flag))
			return n;
		else
			return high;
	}
	int check_leaf(RangeTree_Node * leaf, int xmin, int xmax, int ymin, int ymax, bool xmin_flag = NEG_INFINITY, bool xmax_flag = POS_INFINITY, bool ymin_flag = NEG_INFINITY, bool ymax_flag = POS_INFINITY) {
		if ((!pts_x_smaller(*leaf->split, xmin, xmin_flag)) && pts_x_smaller(*leaf->split, xmax, xmax_flag) &&
			(!pts_y_smaller(*leaf->split, ymin, ymin_flag)) && pts_y_smaller(*leaf->split, ymax, ymax_flag))
			return 1;
		else
			return 0;
	}
	void query_y_window(RangeTree_Node * split_node, int ymin, int ymax, int &ymin_index, int &ymax_index, bool ymin_flag = NEG_INFINITY, bool ymax_flag = POS_INFINITY) {
		int pt_size = split_node->points_and_link.size();
		ymin_index = binary_search_y_first_gt(split_node->points_and_link, ymin, ymin_flag, pt_size);
		ymax_index = binary_search_y_first_gt(split_node->points_and_link, ymax, ymax_flag, pt_size);
	}
	void cascade_to_child(RangeTree_Node * root, bool direction, int &ymin_index, int &ymax_index) {
		// that child cannot be nullptr

		int head_y_size = root->points_and_link.size();
		if (direction == LEFT) {
			RangeTree_Node *p = root->left_child;
			int p_y_size = p->points_and_link.size();

			if (ymin_index >= head_y_size) ymin_index = p_y_size;
			else ymin_index = root->points_and_link[ymin_index].left_child_index;

			if (ymax_index >= head_y_size) ymax_index = p_y_size;
			else ymax_index = root->points_and_link[ymax_index].left_child_index;
		}
		else {
			// direction == RIGHT
			RangeTree_Node *p = root->right_child;
			int p_y_size = p->points_and_link.size();

			if (ymin_index >= head_y_size) ymin_index = p_y_size;
			else ymin_index = root->points_and_link[ymin_index].right_child_index;

			if (ymax_index >= head_y_size) ymax_index = p_y_size;
			else ymax_index = root->points_and_link[ymax_index].right_child_index;
		}
	}
	int query_window_count(Window &query, bool xmin_tight = true, bool xmax_tight = true, bool ymin_tight = true, bool ymax_tight = true) {
		if (head == nullptr) return 0;
		
		int xmin = query.xmin;
		bool xmin_flag = xmin_tight ? NEG_INFINITY : POS_INFINITY;
		int xmax = query.xmax;
		bool xmax_flag = xmax_tight ? POS_INFINITY : NEG_INFINITY;
		int ymin = query.ymin;
		bool ymin_flag = ymin_tight ? NEG_INFINITY : POS_INFINITY;
		int ymax = query.ymax;
		bool ymax_flag = ymax_tight ? POS_INFINITY : NEG_INFINITY;
		// xmin < xmax
		RangeTree_Node * p = this->head;
		while (!p->is_leaf()) {
			if (!pts_x_smaller(*p->split, xmax, xmax_flag)) {
				p = p->left_child;
				continue;
			}
			if (pts_x_smaller(*p->split, xmin, xmin_flag)) {
				p = p->right_child;
				continue;
			}
			
			break;
		}
		RangeTree_Node * split_node = p;

		int result = 0;

		if (split_node->is_leaf()) {
			result += check_leaf(split_node, xmin, xmax, ymin, ymax,
				xmin_flag, xmax_flag, ymin_flag, ymax_flag);
		}
		else {
			int split_ymin_index, split_ymax_index;
			query_y_window(split_node, ymin, ymax, split_ymin_index, split_ymax_index,
				ymin_flag, ymax_flag);

			int v_ymin_index = split_ymin_index;
			int v_ymax_index = split_ymax_index;

			// traverse xmin
			RangeTree_Node * v = split_node->left_child;
			cascade_to_child(split_node, LEFT, v_ymin_index, v_ymax_index);
			while (!v->is_leaf()) {
				if (!pts_x_smaller(*v->split, xmin, xmin_flag)) {
					// report right child count;
					int tmp_ymin_index = v_ymin_index;
					int tmp_ymax_index = v_ymax_index;
					cascade_to_child(v, RIGHT, tmp_ymin_index, tmp_ymax_index);
					result += tmp_ymax_index - tmp_ymin_index;

					// go to left child
					cascade_to_child(v, LEFT, v_ymin_index, v_ymax_index);
					v = v->left_child;
				}
				else {
					// go to right child
					cascade_to_child(v, RIGHT, v_ymin_index, v_ymax_index);
					v = v->right_child;
				}
			}
			result += check_leaf(v, xmin, xmax, ymin, ymax,
				xmin_flag, xmax_flag, ymin_flag, ymax_flag);

			// traverse xmax
			v_ymin_index = split_ymin_index;
			v_ymax_index = split_ymax_index;
			v = split_node->right_child;
			cascade_to_child(split_node, RIGHT, v_ymin_index, v_ymax_index);
			while (!v->is_leaf()) {
				if (pts_x_smaller(*v->split, xmax, xmax_flag)) {
					// report left child count;
					int tmp_ymin_index = v_ymin_index;
					int tmp_ymax_index = v_ymax_index;
					cascade_to_child(v, LEFT, tmp_ymin_index, tmp_ymax_index);
					result += tmp_ymax_index - tmp_ymin_index;

					// go to right child
					cascade_to_child(v, RIGHT, v_ymin_index, v_ymax_index);
					v = v->right_child;
				}
				else {
					// go to left child
					cascade_to_child(v, LEFT, v_ymin_index, v_ymax_index);
					v = v->left_child;
				}
			}
			result += check_leaf(v, xmin, xmax, ymin, ymax,
				xmin_flag, xmax_flag, ymin_flag, ymax_flag);
		}

		return result;
	}
};

class IntervalTree_Node {
public:
	int x_mid;
	IntervalTree_Node * left_child = nullptr;
	IntervalTree_Node * right_child = nullptr;
	RangeTree left_node_range_tree;
	RangeTree right_node_range_tree;
public:
	Point * left_pts = nullptr;
	int left_pts_count = 0;
	Point * right_pts = nullptr;
	int right_pts_count = 0;
public:
	IntervalTree_Node() {}
	~IntervalTree_Node() {}
};

class Point_On_Segment {
public:
	int x;
	int y;
	int segment_id = -1;
public:
	Point_On_Segment() {}
	Point_On_Segment(int mx, int my, int mid) { x = mx; y = my; segment_id = mid; }
	~Point_On_Segment() {}
	const Point_On_Segment& operator =(const Point_On_Segment& b) {
		x = b.x;
		y = b.y;
		segment_id = b.segment_id;
		return b;
	}
};

bool pts_on_segment_smaller(Point_On_Segment &a, Point_On_Segment &b) {
	if (a.x < b.x) return true;
	else if (a.x > b.x) return false;
	else return a.y < b.y;
}

class IntervalTree {
public:
	IntervalTree_Node * head = nullptr;
	int n;
public:
	// temporary
	int *segment_end_point_appears_count;
	Window left_window;
	Window right_window;
public:
	IntervalTree() {}
	~IntervalTree() {}
public:
	//build
	IntervalTree_Node * build_interval_tree(vector<Point_On_Segment> &pts) {
		if (pts.size() == 0) return nullptr;

		int points_count = pts.size();
		int segments_count = pts.size() / 2;
		int mid = segments_count; // 0,1,2,3...mid include mid points
		int x_mid = pts[mid].x;

		IntervalTree_Node * root = new IntervalTree_Node;
		root->x_mid = x_mid;

		int left_segments_count = 0;
		int mid_segments_count = 0;
		int right_segments_count = 0;
		int i = 0;
		
		// count I_left & I_mid
		for (i = 0; i < points_count; i++) {
			if (pts[i].x <= x_mid) {
				if (this->segment_end_point_appears_count[pts[i].segment_id] == 0) {
					this->segment_end_point_appears_count[pts[i].segment_id] = 1;
					mid_segments_count++;
				}
				else if (this->segment_end_point_appears_count[pts[i].segment_id] == 1) {
					//this->segment_end_point_appears_count[pts[i].segment_id] = 2;
					if (pts[i].x < x_mid) {
						// I_left
						mid_segments_count--;
						left_segments_count++;
						this->segment_end_point_appears_count[pts[i].segment_id] = 0;
					}
					else {
						// I_mid
						this->segment_end_point_appears_count[pts[i].segment_id] = 1;
					}
				}
			}
			else break;
		}

		right_segments_count = segments_count - left_segments_count - mid_segments_count;
		vector<Point_On_Segment> left_segments_points(left_segments_count * 2);
		int left_segments_points_i = 0;
		vector<Point_On_Segment> right_segments_points(right_segments_count * 2);
		int right_segments_points_i = 0;

		root->left_pts = new Point[mid_segments_count];
		root->left_pts_count = 0;
		root->right_pts = new Point[mid_segments_count];
		root->right_pts_count = 0;

		for (i = 0; i < points_count; i++) {
			if (pts[i].x <= x_mid) {
				// I_left left&right end points
				if (this->segment_end_point_appears_count[pts[i].segment_id] == 0) {
					left_segments_points[left_segments_points_i] = pts[i];
					left_segments_points_i++;
				}
				// I_mid left end points
				else if (this->segment_end_point_appears_count[pts[i].segment_id] == 1) {
					root->left_pts[root->left_pts_count] = Point(pts[i].x, pts[i].y);
					root->left_pts_count++;
					this->segment_end_point_appears_count[pts[i].segment_id] = 2;
				}
				// I_mid right end points
				else if (this->segment_end_point_appears_count[pts[i].segment_id] == 2) {
					root->right_pts[root->right_pts_count] = Point(pts[i].x, pts[i].y);
					root->right_pts_count++;
					this->segment_end_point_appears_count[pts[i].segment_id] = 3;
				}
			}
			else break;
		}

		// pts[i].x > x_mid
		for (; i < points_count; i++) {
			if (this->segment_end_point_appears_count[pts[i].segment_id] == 0) {
				// I_right left&right end points
				right_segments_points[right_segments_points_i] = pts[i];
				right_segments_points_i++;
			}
			else {
				if (this->segment_end_point_appears_count[pts[i].segment_id] != 2) abort();
				// this->segment_end_point_appears_count[pts[i].segment_id] == 2 
				// I_mid right end points
				root->right_pts[root->right_pts_count] = Point(pts[i].x, pts[i].y);
				root->right_pts_count++;
				this->segment_end_point_appears_count[pts[i].segment_id] = 3;
			}
		}
		
		// build L-tree & R-tree
		root->left_node_range_tree.build_range_tree_from_pts(root->left_pts, root->left_pts_count, false);
		root->right_node_range_tree.build_range_tree_from_pts(root->right_pts, root->right_pts_count, false);

		IntervalTree_Node * root_left_child = build_interval_tree(left_segments_points);
		IntervalTree_Node * root_right_child = build_interval_tree(right_segments_points);
		root->left_child = root_left_child;
		root->right_child = root_right_child;

		return root;
	}
	void build_interval_tree_from_pts(Point_On_Segment * pts, int mn) {
		this->n = mn;
		sort(pts, pts + (mn * 2), pts_on_segment_smaller);

		vector<Point_On_Segment> sorted_points;
		for (int i = 0; i < mn * 2; i++) {
			sorted_points.push_back(pts[i]);
		}

		this->segment_end_point_appears_count = new int[mn];
		for (int i = 0; i < mn; i++) {
			this->segment_end_point_appears_count[i] = 0;
		}
		this->head = build_interval_tree(sorted_points);

		//delete[] this->segment_end_point_appears_count;
		//this->segment_end_point_appears_count = nullptr;
	}
public:
	//query
	int stabbing_query_count(IntervalTree_Node * root, int qx, int ymin, int ymax) {
		if (root == nullptr) return 0;

		// left end point can be on the segment
		int result = 0;
		if (qx < root->x_mid) {
			result += root->left_node_range_tree.query_window_count(this->left_window, true, true, true, true);
			result += stabbing_query_count(root->left_child, qx, ymin, ymax);
		}
		else if (qx > root->x_mid) {
			result += root->right_node_range_tree.query_window_count(this->right_window, true, true, true, true);
			result += stabbing_query_count(root->right_child, qx, ymin, ymax);
		}
		else {
			// qx == root->x_mid
			result += root->left_node_range_tree.query_window_count(this->left_window, true, true, true, true);
		}
		return result;
	}
	int stabbing_query_segment_count(int qx, int ymin, int ymax) {
		if (head == nullptr) return 0;

		this->left_window.xmin = INT_MIN;
		this->left_window.xmax = qx;
		this->left_window.ymin = ymin;
		this->left_window.ymax = ymax;

		this->right_window.xmin = qx;
		this->right_window.xmax = INT_MAX;
		this->right_window.ymin = ymin;
		this->right_window.ymax = ymax;

		return stabbing_query_count(this->head, qx, ymin, ymax);
	}
};

class Search_Structure {
public:
	IntervalTree horizontal_interval_tree;
	RangeTree horizontal_left_end_points_range_tree;
	IntervalTree vertical_interval_tree;
	RangeTree vertical_left_end_points_range_tree;
public:
	Search_Structure() {}
	~Search_Structure() {}
	void build(Point_On_Segment * horizontal_segment_points, Point * horizontal_segment_left_end_points, int horizontal_segment_count,
		Point_On_Segment * vertical_segment_points, Point * vertical_segment_left_end_points, int vertical_segment_count ) {
		if (horizontal_segment_count > 0) {
			this->horizontal_interval_tree.build_interval_tree_from_pts(horizontal_segment_points, horizontal_segment_count);
			this->horizontal_left_end_points_range_tree.build_range_tree_from_pts(horizontal_segment_left_end_points, horizontal_segment_count);
		}
		if (vertical_segment_count > 0) {
			this->vertical_interval_tree.build_interval_tree_from_pts(vertical_segment_points, vertical_segment_count);
			this->vertical_left_end_points_range_tree.build_range_tree_from_pts(vertical_segment_left_end_points, vertical_segment_count);
		}
	}
	int query_window_count(Window &window) {
		int a = horizontal_left_end_points_range_tree.query_window_count(window, false, true, true, true);
		int b = horizontal_interval_tree.stabbing_query_segment_count(window.xmin, window.ymin, window.ymax);
		Window rotate_window;
		rotate_window.xmin = window.ymin;
		rotate_window.xmax = window.ymax;
		rotate_window.ymin = window.xmin;
		rotate_window.ymax = window.xmax;
		int c = vertical_left_end_points_range_tree.query_window_count(rotate_window, false, true, true, true);
		int d = vertical_interval_tree.stabbing_query_segment_count(rotate_window.xmin, rotate_window.ymin, rotate_window.ymax);
		int result = a + b + c + d;
		return result;
	}
};

int main() {
	int n, m;
	scanf("%d %d", &n, &m);
	Point * pts = new Point[n];
	int x1, x2, y1, y2;

	Point_On_Segment * horizontal_segment_points = new Point_On_Segment[n * 2];
	int horizontal_segment_count = 0;
	Point * horizontal_segment_left_end_points = new Point[n];

	Point_On_Segment * vertical_segment_points = new Point_On_Segment[n * 2];
	int vertical_segment_count = 0;
	Point * vertical_segment_left_end_points = new Point[n];
	for (int i = 0; i < n; i++) {
		scanf("%d %d %d %d", &x1, &y1, &x2, &y2);
		if (y1 == y2) {
			// horizontal segment
			if (x1 > x2) swap(x1, x2);
			horizontal_segment_points[horizontal_segment_count * 2].x = x1;
			horizontal_segment_points[horizontal_segment_count * 2].y = y1;
			horizontal_segment_points[horizontal_segment_count * 2].segment_id = horizontal_segment_count;

			horizontal_segment_left_end_points[horizontal_segment_count].x = x1;
			horizontal_segment_left_end_points[horizontal_segment_count].y = y1;

			horizontal_segment_points[horizontal_segment_count * 2 + 1].x = x2;
			horizontal_segment_points[horizontal_segment_count * 2 + 1].y = y1;
			horizontal_segment_points[horizontal_segment_count * 2 + 1].segment_id = horizontal_segment_count;

			horizontal_segment_count++;
		}
		else {
			// x1 == x2
			// vertical segment
			if (y1 > y2) swap(y1, y2);
			vertical_segment_points[vertical_segment_count * 2].x = y1;
			vertical_segment_points[vertical_segment_count * 2].y = x1;
			vertical_segment_points[vertical_segment_count * 2].segment_id = vertical_segment_count;

			vertical_segment_left_end_points[vertical_segment_count].x = y1;
			vertical_segment_left_end_points[vertical_segment_count].y = x1;

			vertical_segment_points[vertical_segment_count * 2 + 1].x = y2;
			vertical_segment_points[vertical_segment_count * 2 + 1].y = x1;
			vertical_segment_points[vertical_segment_count * 2 + 1].segment_id = vertical_segment_count;

			vertical_segment_count++;
		}
	}

	Search_Structure search_structure;
	search_structure.build(horizontal_segment_points, horizontal_segment_left_end_points, horizontal_segment_count,
		vertical_segment_points, vertical_segment_left_end_points, vertical_segment_count);

	Window * windows = new Window[m];
	int * results = new int[m];
	for (int i = 0; i < m; i++) {
		scanf("%d %d %d %d", &windows[i].xmin, &windows[i].ymin, &windows[i].xmax, &windows[i].ymax);
		// query
		results[i] = search_structure.query_window_count(windows[i]);
	}
	for (int i = 0; i < m; i++) {
		printf("%d\n", results[i]);
	}

	return 0;
}
