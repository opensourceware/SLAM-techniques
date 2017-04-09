#include <bits/stdc++.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef std::vector<float> float_vec;
typedef std::pair<float, float> coord_point;

using namespace std;

struct node
{
	float point[2];
	int axis, dim;
	node* left, *right;
};

struct node* newNode(float *point, const int k)
{
	struct node* temp = new node;
	temp->point[0] = point[0];
	temp->point[1] = point[1];
	temp->left=temp->right=NULL;
	temp->dim=k;
};

struct node* insertNode(struct node* root, struct node* n, int axis) {
	if (root==NULL) {
		n->axis = (axis+1)%(n->dim);
		return n;
	}
	if (n->point[root->axis] <= root->point[root->axis]){
		root->left = insertNode(root->left, n, root->axis);
	}
	else {
		root->right = insertNode(root->right, n, root->axis);
	}
	return root;
}

struct node* constructKDTree(std::vector<coord_point> p, const int k) {
	if (p.empty()) {
		return NULL;
	}
	struct node* root;
	int i=0;
	for (coord_point coord: p) {
	//for (int i=0; i<p.size(); i++) {
		float point[] = {coord.first, coord.second};
	//	float point[] = {p[i].first, p[i].second};
		if (i==0) {
			i++;
			root = newNode(point, k);
			root->axis = 0;
		}
		else {
			struct node* n = newNode(point, k);
			root = insertNode(root, n, root->axis);
			//cout << root->point[0] << " " << root->point[1] << "\n";
		}
	}
	return root;
}

void updateNearestNeighbour(float query[], struct node* root, float NNpoint[], float* distance) {
	//L2 Norm
	int n = sizeof(query)/sizeof(int);
	int sum = 0;
	for (int i=0; i<n; i++) {
		sum += pow(query[i]-root->point[i], 2);
	}
	int distance_ = sqrt(sum);
	if (distance_	 < *distance) {
		NNpoint[0] = root->point[0];
		NNpoint[1] = root->point[1];
		*distance = distance_;
	}
}

void NearestNeighborSearch(float query[], struct node* root, float NNpoint[], float* distance) {
	if (root->left==NULL && root->right==NULL) {
		updateNearestNeighbour(query, root, NNpoint, distance);
	}
	else {
		string s;
		if (query[root->axis] <= root->point[root->axis])
			s = "search_left_first";
		else
			s = "search_right_first";
		if (s=="search_left_first") {
			//Search left
			if (query[root->axis]-*distance <= root->point[root->axis]) {
				if (root->left == NULL)
					updateNearestNeighbour(query, root, NNpoint, distance);
				else
					NearestNeighborSearch(query, root->left, NNpoint, distance);
			}
			//Search right
			if (query[root->axis]+*distance > root->point[root->axis]) {
				if (root->right == NULL)
					updateNearestNeighbour(query, root, NNpoint, distance);
				else
					NearestNeighborSearch(query, root->right, NNpoint, distance);
			}
		}
		else //s=="search_right_first" 
		{
			if (query[root->axis]+*distance > root->point[root->axis]) {
				if (root->right == NULL)
					updateNearestNeighbour(query, root, NNpoint, distance);
				else
					NearestNeighborSearch(query, root->right, NNpoint, distance);
			}
			if (query[root->axis]-*distance <= root->point[root->axis]) {
				if (root->left == NULL)
					updateNearestNeighbour(query, root, NNpoint, distance);
				else
					NearestNeighborSearch(query, root->left, NNpoint, distance);
			}
		}
	}
	return;
}

void printTree(struct node* root, int indent) {
	if (root!=NULL) {
		cout << root->point[0] << ' ' << root->point[1] << "\n";
		if (root->left) printTree(root->left, indent+4);
		if (root->right) printTree(root->right, indent+4);
		if (indent) {
			std::cout << std::setw(indent) << ' ';
		}
	}
}


/*
int main() {
	std::vector<coord_point> p(4);
	p[0] = make_pair(5.0, 5.0);
	p[1] = make_pair(0.0, 8.0);
	p[2] = make_pair(7.0, 0.0);
	p[3] = make_pair(7.0, 8.0);
	const int k = 2;

	struct node* root = constructKDTree(p, k);
	//printTree(root, 0);
	float query[2] = {8.0, 8.0};
	float neighbor[2] = {root->point[0], root->point[1]};
	//float* distance;
	float dist = 4.24;
	float* distance = &dist;
	NearestNeighborSearch(query, root, neighbor, distance);
	cout << neighbor[0] << ' ' << neighbor[1];
	cout << *distance;
}

*/