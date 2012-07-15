///////////////////////////////////////////////////////////////////////////
//
// NAME
//  FeatureAlign.h -- image registration using feature matching
//
// SEE ALSO
//  FeatureAlign.h      longer description
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "FeatureAlign.h"
#include "P3Math.h"
#include <math.h>
#include <set>
#include <time.h>
#include <iostream>


#define LargeNum 1000000

/******************* TO DO *********************
 * alignPair:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *		m: motion model
 *		f: focal length
 *		width: image width
 *		height: image height
 *		nRANSAC: number of RANSAC iterations
 *		RANSACthresh: RANSAC distance threshold
 *		M: transformation matrix (output)
 *	OUTPUT:
 *		repeat for nRANSAC iterations:
 *			choose a minimal set of feature matches
 *			estimate the transformation implied by these matches
 *			count the number of inliers
 *		for the transformation with the maximum number of inliers,
 *		compute the least squares motion estimate using the inliers,
 *		and store it in M
 */

vector<int> RandGen(const vector<FeatureMatch> &matches)
{
	srand((unsigned)time(NULL));
	vector<int> vrandNums;
	for(int i = 0; i < 1e4; i++)
	{
		vrandNums.push_back(rand() % matches.size());
	}
	return vrandNums;
}

int RandSelectInliers(const vector<int>& vrandNums, int offPos,
					   const vector<FeatureMatch>& matches,
					   vector<int>& inliers, 
					   set<int>& visitedPair)
{
	inliers.clear();
	int fpos1, fpos2;
	for(int i = offPos; i < (int)vrandNums.size(); i+=2)
	{
		fpos1 = vrandNums[i];
		fpos2 = vrandNums[i+1];
		
		set<int>::iterator iter = visitedPair.find(fpos1 * LargeNum + fpos2);
		if((iter != visitedPair.end()) || (fpos1 == fpos2) || 
			(matches[fpos1].id < 0) || (matches[fpos2].id < 0))
		{
			continue;
		}
		else
		{
			inliers.push_back(fpos1);
			inliers.push_back(fpos2);
			visitedPair.insert(fpos1 * LargeNum + fpos2);
			break;
		}
	} 

	return i+2;
}

vector<int> alignPair(const FeatureSet &f1, const FeatureSet &f2,
			  const vector<FeatureMatch> &matches, MotionModel m,
			  float f, const pair<int,int>& WH1, const pair<int,int>& WH2,
			  int nRANSAC, double RANSACthresh, CTransform3x3& M)
{
	// BEGIN TODO
	// write this entire method
	vector< vector<int> > vInliers;
	set<int> visitedPair;
	vector<int> vRandNums = RandGen(matches);
	int offPos = 0;

	for(int i = 0; i < nRANSAC; i++)
	{
		vector<int> inliers;
		offPos = RandSelectInliers(vRandNums, offPos, matches, inliers, visitedPair);
		if(inliers.size() < 2)
			break;

		CTransform3x3 Mtmp;
		leastSquaresFit(f1, f2, matches, m, f, WH1, WH2, inliers, Mtmp);

		countInliers(f1,f2,matches,m,f,WH1,WH2,Mtmp,RANSACthresh,inliers);

		vInliers.push_back(inliers);
	}

	int tmax = -1, maxPos;
	for (size_t j = 0; j < vInliers.size(); j++)
	{
		if(tmax < (int)vInliers[j].size())
		{
			tmax = vInliers[j].size();
			maxPos = j;
		}
	}

	if(vInliers[maxPos].size() < 2)
		return vInliers[maxPos];

	leastSquaresFit(f1, f2, matches, m,f,WH1,WH2,vInliers[maxPos],M);
	// END TODO

	return vInliers[maxPos];
}

void PrintVector(const CVector3 &v)
{
	cout<<endl;
	cout<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<endl;
}

/******************* TO DO *********************
 * countInliers:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *		m: motion model
 *		f: focal length
 *		width: image width
 *		height: image height
 *		M: transformation matrix
 *		RANSACthresh: RANSAC distance threshold
 *		inliers: inlier feature IDs
 *	OUTPUT:
 *		transform the features in f1 by M
 *
 *		count the number of features in f1 for which the transformed
 *		feature is within Euclidean distance RANSACthresh of its match
 *		in f2
 *
 *		store these features IDs in inliers
 *
 *		this method should be similar to evaluateMatch from project 1,
 *		except you are comparing each distance to a threshold instead
 *		of averaging them
 */
int countInliers(const FeatureSet &f1, const FeatureSet &f2,
				 const vector<FeatureMatch> &matches, MotionModel m,
				 float f, const pair<int,int>& WH1, const pair<int,int>& WH2,
				 CTransform3x3 M, double RANSACthresh, vector<int> &inliers)
{
	inliers.clear();
	int count = 0;

	for (unsigned int i=0; i<f1.size(); i++) {
		// BEGIN TODO
		// determine if the ith feature in f1, when transformed by M,
		// is within RANSACthresh of its match in f2 (if one exists)
		//
		// if so, increment count and append i to inliers
		if(matches[i].id < 0)
			continue;
		
		CVector3 p,q;
		//p[0] = f1[i].x;	p[1] = f1[i].y;	p[2] = f;
		int width = WH1.first,	height = WH1.second;
		p[0] = f1[i].x - width / 2.0;	p[1] = f1[i].y - height / 2.0;	p[2] = f;
		q = M * p;
		//PrintVector(q);
		
		//double xNew = q[0];
		//double yNew = q[1];
		width = WH2.first;	height = WH2.second;
		double xNew = q[0] + width / 2.0;
		double yNew = q[1] + height / 2.0;
		int f2pos = matches[i].id - 1;
		double dist = sqrt(pow(xNew-f2[f2pos].x,2) + pow(yNew-f2[f2pos].y,2));
		if(dist < RANSACthresh)
		{
			inliers.push_back(i);
			count++; 
		}
		// END TODO
	}

	return count;
}


void NormRotationMat(CTransform3x3 &M)
{
	CVector3 p,q,r,vtmp;
	p[0]=M[0][0];	p[1]=M[1][0];	p[2]=M[2][0];
	q[0]=M[0][1];	q[1]=M[1][1];	q[2]=M[2][1];
	r[0]=M[0][2];	r[1]=M[1][2];	r[2]=M[2][2];
	p.Normalize();
	q.Normalize();
	r.Normalize();
	vtmp = p.cross(q);
	double dist = sqrt(pow(vtmp[0]-r[0],2)+pow(vtmp[1]-r[1],2)+pow(vtmp[2]-r[2],2));
	if(dist > 1.0e-6)
		r = vtmp;
	M[0][0] = p[0];	M[0][1] = q[0];	M[0][2] = r[0];
	M[1][0] = p[1];	M[1][1] = q[1];	M[1][2] = r[1];
	M[2][0] = p[2];	M[2][1] = q[2];	M[2][2]	= r[2];
}

void PrintMat(const CTransform3x3 &M)
{
	cout<<endl;
	cout<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[0][2]<<endl;
	cout<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<M[1][2]<<endl;
	cout<<M[2][0]<<"\t"<<M[2][1]<<"\t"<<M[2][2]<<endl;
}

/******************* TO DO *********************
 * leastSquaresFit:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *		m: motion model
 *		f: focal length
 *		width: image width
 *		height: image height
 *		inliers: inlier feature IDs
 *		M: transformation matrix (output)
 *	OUTPUT:
 *		compute the transformation from f1 to f2 using only the inliers
 *		and return it in M
 */
int leastSquaresFit(const FeatureSet &f1, const FeatureSet &f2,
					const vector<FeatureMatch> &matches, MotionModel m,
					float f, const pair<int,int>& WH1, const pair<int,int>& WH2,
					const vector<int> &inliers, CTransform3x3& M)
{
	// BEGIN TODO
	// write this entire method
	
	CTransform3x3 mat(0), tempMat(0);
	CVector3 v1, v2;
	
	for(size_t i = 0; i < inliers.size(); i++)
	{
		int f1pos = inliers[i];
		int f2pos = matches[f1pos].id - 1;
		//v1[0] = f1[f1pos].x;	v1[1] = f1[f1pos].y;	v1[2] = f;
		//v2[0] = f2[f2pos].x;	v2[1] = f2[f2pos].y;	v2[2] = f;
		int width = WH1.first, height = WH1.second;
		v1[0] = f1[f1pos].x-width/2.0;	v1[1] = f1[f1pos].y-height/2.0;	v1[2] = f;
		width = WH2.first;	height = WH2.second;
		v2[0] = f2[f2pos].x-width/2.0;	v2[1] = f2[f2pos].y-height/2.0;	v2[2] = f;
		tempMat = v1 * v2;
		mat = mat + tempMat;
	}

	CTransform3x3 u, s, v;
	svd(mat, u, s, v);

	M = v * u.Transpose();

	//PrintMat(M);
	//NormRotationMat(M);
	//PrintMat(M);
	// END TODO

	return 0;
}
