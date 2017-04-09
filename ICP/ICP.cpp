#include <bits/stdc++.h>
//#include <fstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include </usr/local/include/eigen3/Eigen/Dense>
#include "/home/manpreet/Desktop/UMass/Robotics/assignment 3/GraphSLAM/kd-trees.cpp"
#include <iostream>
#include <fstream>


using namespace std;
//using namespace stdext;
using namespace Eigen;

typedef std::pair<float, float> coord_point;
typedef std::map<float*, float*> point_map;
typedef std::vector<float> float_vec;
typedef std::vector<std::pair<coord_point, coord_point>> floatp_vec;
typedef std::pair<float*, float*> coord_pointp;
typedef std::pair<coord_point, coord_point> coord_pair;

#define PI 3.14159265

float convert_to_rad(float angle) {
	return angle*PI/180;
}

std::vector<coord_point> convert2coord(float_vec scan, float_vec pose, float_vec laser_theta) {
	float x_pose = pose[0];
	float y_pose = pose[1];
	float theta_pose = pose[2];
	std::vector<coord_point> coord(scan.size());
	float cos_theta_pose;
	float sin_theta_pose;
	int i=0;
	for (float f: scan) {
		float cos_theta_pose = cos(convert_to_rad(theta_pose+laser_theta[i]));
		float sin_theta_pose = sin(convert_to_rad(theta_pose+laser_theta[i]));
		coord[i++] = std::make_pair(x_pose+f*cos_theta_pose, y_pose+f*sin_theta_pose);
	}
	return coord;
}

floatp_vec icp_map(std::vector<coord_point> reference_points, std::vector<coord_point> source_points)	{
	struct node* kdtree_root = constructKDTree(reference_points, 2);
	float query[2];
	float* distance;
	float dist = 50;
	distance = &dist;
	//point_map mapping;
	floatp_vec mapping(180);
	int i = 0;
	for (coord_point source_point: source_points) {
		query[0] = source_point.first, 
		query[1] = source_point.second;
		float neighbor[2] = {kdtree_root->point[0], kdtree_root->point[1]};
		dist = sqrt(pow((query[0]-neighbor[0]), 2)+pow((query[1]-neighbor[1]), 2));
		NearestNeighborSearch(query, kdtree_root, neighbor, distance);
		mapping[i++] = make_pair(make_pair(query[0], query[1]), make_pair(neighbor[0], neighbor[1]));
	}
	return mapping;
}

MatrixXf lsqmean(std::vector<coord_point> &source, std::vector<coord_point> reference, std::vector<coord_point> originalReference, std::string filename) 
{
	MatrixXf Rot = MatrixXf::Random(3, 2);
	float err = 80000;
	MatrixXf OriginalSrc = MatrixXf::Random(source.size(), 2);
	MatrixXf OriginalRef = MatrixXf::Random(source.size(), 2);
	MatrixXf Src = MatrixXf::Random(source.size(), 3);
	MatrixXf Src_ = MatrixXf::Random(source.size(), 2);
	MatrixXf Ref = MatrixXf::Random(source.size(), 2);
	for (int i=0; i<source.size(); i++) {
		OriginalSrc.row(i) << source[i].first, source[i].second;
		OriginalRef.row(i) << originalReference[i].first, originalReference[i].second;
		Src.row(i) << source[i].first, source[i].second, 1.0;
		Ref.row(i) << reference[i].first, reference[i].second;
	}
	while (err > 700) {
		Rot = Src.jacobiSvd(ComputeThinU | ComputeThinV).solve(Ref);
		MatrixXf AbsError = (Src*Rot-Ref).cwiseAbs();
		AbsError = AbsError.cwiseProduct(AbsError);
		err = AbsError.sum();
		cout<< err << "\n";
		Src_ = Src*Rot;
		for (int i=0; i<source.size(); i++) {
			Src.row(i) << Src_(i, 0), Src_(i, 1), 1.0;
		}
		//break;
	}
	for (int i=0; i<source.size(); i++) {
		source[i] = make_pair(Src_(i, 0), Src_(i, 1));
		//cout << source[i].first << " " << source[i].second << " " << Src_(i, 0) << " " << Src_(i, 1) << "\n";
	}
	std::ofstream file(filename);
	if (file.is_open())
	{
		file << "Here is the matrix Src:\n" << OriginalSrc << '\n';
		file << "Here is the matrix Ref:\n" << Ref << '\n';
		file << "Here is the matrix Src_corrected:\n" << Src_ << '\n';
		file << "Here is the matrix OriginalRef:\n" << OriginalRef << '\n';
		file << "Here is the Rotation Matrix:\n" << Rot << "\n";
	}
//	}
	return Rot;
}

int main()
{
	std::ifstream input("intel.script");
	char* pch;
	int flag = 1;
	string pos = "POS";
	string laser = "LASER-RANGE";
	std::vector<float> reference_scan(180);
	std::vector<float> source_scan(180);
	std::vector<coord_point> reference_points(180);
	std::vector<coord_point> source_points(180);
	std::vector<coord_point> mapped_reference_points(180);
	std::vector<float> o(3);
	std::vector<float> prev_state(3);
	std::vector<float> state(3);
	state[0] = 0.0;
	state[1] = 0.0;
	state[2] = 0.0;
	int line_num = 0;
	for (std::string line; std::getline(input, line);) {
		if (line_num<500) {
			line_num++;
			continue;
		}
		else 
			line_num++;
		if (line_num==540)
			break;
		char *str = new char[line.size()+1];
		std::copy(line.begin(), line.end(), str);
		str[line.size()] = '\0';
		pch = strtok(str, " :");
		if (pch==pos) {
			int i=0;
			while(pch!=NULL)
			{
				if (i<=2) {
					i++;
				}
				else {
					double num = ::atof(pch);
					o[(i++)-3] = num;
				}
				pch = std::strtok(NULL, " :");
			}
			for (int j=0; j<=3; j++) {
				prev_state[j] = state[j];
				state[j] = o[j];
			}
		}
		else if (pch==laser) {
			int i=0;
			float minangle;
			float num_scans;
			float maxangle;
			while(pch!=NULL) {
				if (i==3) {
					minangle = ::atof(pch);
				}
				else if (i==4) {
					num_scans = ::atof(pch);
				}
				else if (i==5) {
					maxangle = ::atof(pch);
				}
				else if (i>5) {
					float value = ::atof(pch);
					source_scan[i-6] = value;
				}
				++i;
				pch = std::strtok(NULL, " :");
			}	
			float increment = (maxangle-minangle)/num_scans;
			float_vec laser_theta(num_scans);
			for (int i=0; i<laser_theta.size(); i++) {
				laser_theta[i] = minangle+i*increment-90	;
			}
			if (flag) {
				for (int j=0; j<180; j++) {
					reference_scan[j] = source_scan[j];
					reference_points = convert2coord(reference_scan, prev_state, laser_theta);
				}
				flag = 0;
			}
			else {
				source_points = convert2coord(source_scan, state, laser_theta);
				floatp_vec icp_mapping = icp_map(reference_points, source_points);
				i = 0;
				for (coord_pair m: icp_mapping) {
					source_points[i] = m.first;
					mapped_reference_points[i++] = m.second;
				}
				std::string filename = std::to_string(line_num);
				MatrixXf Rot = lsqmean(source_points, mapped_reference_points, reference_points, filename);

				for (int j=0; j<180; j++) {
					reference_scan[i] = source_scan[i];
				}
				//Updates the reference point to the previous timestep
				//If you want to use global registration, comment this next line.
				//reference_points = source_points;
			}
		}
	}
	return 0;
}
