
#include "float.h"
#include "math.h"


#include "ba.hpp"
#include "cali.hpp"
#include "relativepose.hpp"
#include "distortion.hpp"
#include "bundlerio.hpp"

//matrix
#include "matrix/matrix.h"


#ifndef WIN32
#include <ext/hash_set>
#include <ext/hash_map>
#else
#include <hash_set>
#include <hash_map>
#endif

#ifdef WIN32
#define isnan _isnan
#endif


#include <assert.h>
#include <algorithm>

//sfm driver lib
#include "sfm.h"

//imagelib
#include "qsort.h"
#include "util.h"
#include "defines.h"

//matrix
#include "matrix.h"

//triangulate
#include "triangulate.h"

#include "bundlerio.hpp"


static char ply_header1[] = 
"ply\n"
"format ascii 1.0\n"
"element face 0\n"
"property list uchar int vertex_indices\n"
"element vertex %d\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property uchar diffuse_red\n"
"property uchar diffuse_green\n"
"property uchar diffuse_blue\n"
"end_header\n";

/* Write point files to a ply file */
void DumpPointsToPly(char *output_directory, char *filename, 
							  int num_points, int num_cameras, 
							  v3_t *points, v3_t *colors,
							  camera_params_t *cameras 
							  /*bool reflect*/) 
{
	int num_good_pts = 0;

	for (int i = 0; i < num_points; i++) 
	{
		if (Vx(colors[i]) == 0x0 && 
			Vy(colors[i]) == 0x0 && 
			Vz(colors[i]) == 0xff) 
			continue;
		num_good_pts++;
	}

	char ply_out[256];
	sprintf(ply_out, "%s/%s", output_directory, filename);

	FILE *f = fopen(ply_out, "w");

	if (f == NULL) 
	{
		printf("Error opening file %s for writing\n", ply_out);
		return;
	}

	/* Print the ply header */
	fprintf(f, ply_header1, num_good_pts + 2 * num_cameras);

	/* Now triangulate all the correspondences */
	for (int i = 0; i < num_points; i++) 
	{
		if (Vx(colors[i]) == 0x0 && 
			Vy(colors[i]) == 0x0 && 
			Vz(colors[i]) == 0xff) 
			continue;

		/* Output the vertex */
		fprintf(f, "%0.6e %0.6e %0.6e %d %d %d\n", 
			Vx(points[i]), Vy(points[i]), Vz(points[i]),
			// Vx(points[idx]), Vy(points[idx]), Vz(points[idx]),
			// (reflect ? -1 : 1) * Vz(points[i]),
			iround(Vx(colors[i])), 
			iround(Vy(colors[i])), 
			iround(Vz(colors[i])));
	}

	for (int i = 0; i < num_cameras; i++) 
	{
		double c[3];

		double Rinv[9];
		matrix_invert(3, cameras[i].R, Rinv);

		memcpy(c, cameras[i].t, 3 * sizeof(double));

		if ((i % 2) == 0)
			fprintf(f, "%0.6e %0.6e %0.6e 0 255 0\n", c[0], c[1], c[2]);
		// (reflect ? -1 : 1) * c[2]);
		else
			fprintf(f, "%0.6e %0.6e %0.6e 255 0 0\n", c[0], c[1], c[2]);
		// (reflect ? -1 : 1) * c[2]);

		double p_cam[3] = { 0.0, 0.0, -0.05 };
		double p[3];

		// if (!reflect)
		//    p_cam[2] *= -1.0;

		matrix_product(3, 3, 3, 1, Rinv, p_cam, p);

		p[0] += c[0];
		p[1] += c[1];
		p[2] += c[2];

		fprintf(f, "%0.6e %0.6e %0.6e 255 255 0\n",
			p[0], p[1], p[2]); // (reflect ? -1 : 1) * p[2]);
	}
	fclose(f);
}

/* Dump an output file containing information about the current
* state of the world */
void DumpOutputFile(char *output_dir, char *filename, 
							 int num_images, int num_cameras, int num_points,
							 int *added_order, 
							 camera_params_t *cameras, 
							 v3_t *points, v3_t *colors,
							 std::vector<ImageKeyVector> &pt_views
							 /*bool output_radial_distortion*/)
{
	//clock_t start = clock();

	int num_visible_points = 0;

	for (int i = 0; i < num_points; i++) 
	{
		if (pt_views[i].size() > 0)
			num_visible_points++;
	}

	char buf[256];
	sprintf(buf, "%s/%s", output_dir, filename);

	FILE *f = fopen(buf, "w");
	if (f == NULL) 
	{
		printf("Error opening file %s for writing\n", buf);
		return;
	}

	// if (output_radial_distortion) {
	/* Print version number */
	// fprintf(f, "# Bundle file v0.4\n");
	fprintf(f, "# Bundle file v0.3\n");
	// }

	fprintf(f, "%d %d\n", num_images, num_visible_points);

	/* Dump cameras */
	for (int i = 0; i < num_images; i++) 
	{

#if 0
		/* Print the name of the file */
		fprintf(f, "%s %d %d\n", 
			m_image_data[i].m_name, 
			m_image_data[i].GetWidth(), m_image_data[i].GetHeight());
#endif

		int idx = -1;
		for (int j = 0; j < num_cameras; j++) {
			if (added_order[j] == i) {
				idx = j;
				break;
			}
		}

		if (idx == -1) {
			// if (!output_radial_distortion)
			//     fprintf(f, "0\n");
			// else
			fprintf(f, "0 0 0\n");
			fprintf(f, "0 0 0\n0 0 0\n0 0 0\n0 0 0\n");
		} else {
			// if (!output_radial_distortion)
			//     fprintf(f, "%0.10e\n", cameras[idx].f);
			// else
			fprintf(f, "%0.10e %0.10e %0.10e\n", 
				cameras[idx].f, cameras[idx].k[0], cameras[idx].k[1]);

			fprintf(f, "%0.10e %0.10e %0.10e\n", 
				cameras[idx].R[0], 
				cameras[idx].R[1], 
				cameras[idx].R[2]);
			fprintf(f, "%0.10e %0.10e %0.10e\n", 
				cameras[idx].R[3], 
				cameras[idx].R[4], 
				cameras[idx].R[5]);
			fprintf(f, "%0.10e %0.10e %0.10e\n", 
				cameras[idx].R[6], 
				cameras[idx].R[7], 
				cameras[idx].R[8]);

			double t[3];
			matrix_product(3, 3, 3, 1, cameras[idx].R, cameras[idx].t, t);
			matrix_scale(3, 1, t, -1.0, t);
			fprintf(f, "%0.10e %0.10e %0.10e\n", t[0], t[1], t[2]);
		}
		fprintf(f, "\n");
	}

	/* Dump points */
	/*
	for (int i = 0; i < num_points; i++) 
	{
		int num_visible = (int) pt_views[i].size();

		if (num_visible > 0) 
		{
			// Position 
			fprintf(f, "%0.10e %0.10e %0.10e\n", 
				Vx(points[i]), Vy(points[i]), Vz(points[i]));
			// Vx(points[idx]), Vy(points[idx]), Vz(points[idx]));

			// Color
			fprintf(f, "%d %d %d\n", 
				iround(Vx(colors[i])), 
				iround(Vy(colors[i])), 
				iround(Vz(colors[i])));

			int num_visible = (int) pt_views[i].size();
			fprintf(f, "%d", num_visible);
			for (int j = 0; j < num_visible; j++) 
			{
				int img = added_order[pt_views[i][j].first];
				int key = pt_views[i][j].second;

				double x = m_image_data[img].m_keys[key].m_x;
				double y = m_image_data[img].m_keys[key].m_y;
				fprintf(f, " %d %d %0.4f %0.4f", img, key, x, y);
			}

			fprintf(f, "\n");
		}
	}
	*/

#if 0
	/* Finally, dump all outliers */
	ImageKeyVector outliers;
	for (int i = 0; i < num_images; i++) {
		/* Find the index of this camera in the ordering */
		int idx = -1;
		for (int j = 0; j < num_cameras; j++) {
			if (added_order[j] == i) {
				idx = j;
				break;
			}
		}

		if (idx == -1) continue;

		int num_keys = GetNumKeys(i);
		for (int j = 0; j < num_keys; j++) {
			if (GetKey(i,j).m_extra == -2) {
				outliers.push_back(ImageKey(i,j));
			}
		}
	}

	int num_outliers = (int) outliers.size();
	fprintf(f, "%d\n", num_outliers);

	for (int i = 0; i < num_outliers; i++) {
		fprintf(f, "%d %d\n", outliers[i].first, outliers[i].second);
	}
#endif

	fclose(f);

	//clock_t end = clock();

	//printf("[SifterApp::DumpOutputFile] Wrote file in %0.3fs\n",
	//	(double) (end - start) / (double) CLOCKS_PER_SEC);
}


static int compare_doubles(const void *d1, const void *d2)
{
	double a = *(double *) d1;
	double b = *(double *) d2;

	if (a < b) return -1;
	if (a > b) return 1;
	return 0;
}

/* Compute the angle between two rays */
double ComputeRayAngle(v2_t p, v2_t q, 
					   const camera_params_t &cam1, 
					   const camera_params_t &cam2)
{
	double K1[9], K2[9];
	GetIntrinsics(cam1, K1);
	GetIntrinsics(cam2, K2);

	double K1_inv[9], K2_inv[9];
	matrix_invert(3, K1, K1_inv);
	matrix_invert(3, K2, K2_inv);

	double p3[3] = { Vx(p), Vy(p), 1.0 };
	double q3[3] = { Vx(q), Vy(q), 1.0 };

	double p3_norm[3], q3_norm[3];
	matrix_product331(K1_inv, p3, p3_norm);
	matrix_product331(K2_inv, q3, q3_norm);

	v2_t p_norm = v2_new(p3_norm[0] / p3_norm[2], p3_norm[1] / p3_norm[2]);
	v2_t q_norm = v2_new(q3_norm[0] / q3_norm[2], q3_norm[1] / q3_norm[2]);

	double R1_inv[9], R2_inv[9];
	matrix_transpose(3, 3, (double *) cam1.R, R1_inv);
	matrix_transpose(3, 3, (double *) cam2.R, R2_inv);

	double p_w[3], q_w[3];

	double pv[3] = { Vx(p_norm), Vy(p_norm), -1.0 };
	double qv[3] = { Vx(q_norm), Vy(q_norm), -1.0 };

	double Rpv[3], Rqv[3];

	matrix_product331(R1_inv, pv, Rpv);
	matrix_product331(R2_inv, qv, Rqv);

	matrix_sum(3, 1, 3, 1, Rpv, (double *) cam1.t, p_w);
	matrix_sum(3, 1, 3, 1, Rqv, (double *) cam2.t, q_w);

	/* Subtract out the camera center */
	double p_vec[3], q_vec[3];
	matrix_diff(3, 1, 3, 1, p_w, (double *) cam1.t, p_vec);
	matrix_diff(3, 1, 3, 1, q_w, (double *) cam2.t, q_vec);

	/* Compute the angle between the rays */
	double dot;
	matrix_product(1, 3, 3, 1, p_vec, q_vec, &dot);

	double mag = matrix_norm(3, 1, p_vec) * matrix_norm(3, 1, q_vec);

	return acos(CLAMP(dot / mag, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
}





v3_t GeneratePointAtInfinity(vector<CImageDataBase*> imageData,
							 const ImageKeyVector &views, 
							 int *added_order, 
							 camera_params_t *cameras,
							 double &error, 
							 bool explicit_camera_centers)
{
	camera_params_t *cam = NULL;

	int camera_idx = views[0].first;
	int image_idx = added_order[camera_idx];
	int key_idx = views[0].second;
	FeatPoint &key = imageData[image_idx]->GetKeyPoint(key_idx);//GetKey(image_idx, key_idx);

	cam = cameras + camera_idx;

	double p3[3] = { key.cx, key.cy, 1.0 };
	


	double K[9], Kinv[9];
	GetIntrinsics(cameras[camera_idx], K);
	matrix_invert(3, K, Kinv);

	double ray[3];
	matrix_product(3, 3, 3, 1, Kinv, p3, ray);

	//We now have a ray, put it at infinity
	double ray_world[3];
	matrix_transpose_product(3, 3, 3, 1, cam->R, ray, ray_world);

	double pos[3] = { 0.0, 0.0, 0.0 };
	double pt_inf[3] = { 0.0, 0.0, 0.0 };

	if (!explicit_camera_centers) 
	{

	}
	else
	{
		memcpy(pos, cam->t, 3 * sizeof(double));
		double ray_extend[3];
		matrix_scale(3, 1, ray, 100.0, ray_extend);
		matrix_sum(3, 1, 3, 1, pos, ray, pt_inf);
	}

	return v3_new(pt_inf[0], pt_inf[1], pt_inf[2]);
}

// Triangulate a subtrack
v3_t TriangulateNViews(const vector<CImageDataBase*> imageData,
					   const ImageKeyVector &views, 
					   int *added_order, camera_params_t *cameras,
					   double &error, bool explicit_camera_centers)
{
	bool estimate_distortion = true;

	int num_views = (int) views.size();

	v2_t *pv = new v2_t[num_views];
	double *Rs = new double[9 * num_views];
	double *ts = new double[3 * num_views];

	for (int i = 0; i < num_views; i++) 
	{
		camera_params_t *cam = NULL;

		int camera_idx = views[i].first;
		int image_idx = added_order[camera_idx];
		int key_idx = views[i].second;
		
		//Keypoint &key = GetKey(image_idx, key_idx);
		//double p3[3] = { key.m_x, key.m_y, 1.0 };
		FeatPoint &key = imageData[image_idx]->GetKeyPoint(key_idx); //GetKey(image_idx, key_idx);
		double p3[3] = { key.cx, key.cy, 1.0 };

		double K[9], Kinv[9];
		GetIntrinsics(cameras[camera_idx], K);
		matrix_invert(3, K, Kinv);

		double p_n[3];
		matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

		// EDIT!!!
		pv[i] = v2_new(-p_n[0], -p_n[1]);
		pv[i] = UndistortNormalizedPoint(pv[i], cameras[camera_idx]);

		cam = cameras + camera_idx;

		memcpy(Rs + 9 * i, cam->R, 9 * sizeof(double));
		if (!explicit_camera_centers) 
		{
			memcpy(ts + 3 * i, cam->t, 3 * sizeof(double));
		}
		else 
		{
			matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * i);
			matrix_scale(3, 1, ts + 3 * i, -1.0, ts + 3 * i);
		}
	}

	v3_t pt = triangulate_n(num_views, pv, Rs, ts, &error);

	error = 0.0;
	for (int i = 0; i < num_views; i++)
	{
		int camera_idx = views[i].first;
		int image_idx = added_order[camera_idx];
		int key_idx = views[i].second;
		
		//Keypoint &key = GetKey(image_idx, key_idx);
		FeatPoint &key = imageData[image_idx]->GetKeyPoint(key_idx);

		v2_t pr = sfm_project_final(cameras + camera_idx, pt, 
			explicit_camera_centers ? 1 : 0,
			estimate_distortion ? 1 : 0);

	
		double dx = Vx(pr) - key.cx; //key.m_x;
		double dy = Vy(pr) - key.cy; //key.m_y;

		error += dx * dx + dy * dy;
	}

	error = sqrt(error / num_views);

	delete [] pv;
	delete [] Rs;
	delete [] ts;

	return pt;
}



/* Add new points to the bundle adjustment */
int BundleAdjustAddAllNewPoints(int num_points, int num_cameras,
								int *added_order,
								camera_params_t *cameras,
								v3_t *points, v3_t *colors,
								double reference_baseline,
								vector<CImageDataBase*> imageData,
								vector<TrackInfo>& trackSeq,
								vector<ImageKeyVector> &pt_views,
								double max_reprojection_error,
								int min_views)
{
	std::vector<int> track_idxs;
	std::vector<ImageKeyVector> new_tracks;

	// __gnu_cxx::hash_map<int,bool> tracks_seen;
	int num_tracks_total = (int) trackSeq.size();
	int *tracks_seen = new int[num_tracks_total];
	for (int i = 0; i < num_tracks_total; i++) 
	{
		tracks_seen[i] = -1;
	}

	/* Gather up the projections of all the new tracks */
	printf("adding new tracks... \n");
	for (int i = 0; i < num_cameras; i++) 
	{
		int image_idx1 = added_order[i];

		int num_keys = imageData[image_idx1]->GetFeatNumber();//GetNumKeys(image_idx1);

		for (int j = 0; j < num_keys; j++) 
		{
			//Keypoint &key = GetKey(image_idx1, j);
			FeatPoint &key = imageData[image_idx1]->GetKeyPoint(j); 

			if (key.track == -1)
				continue;  //Key belongs to no track 

			if (key.extra != -1)
				continue;  // Key is outlier or has already been added 

			int track_idx = key.track;

			// Check if this track is already associated with a point
			if (trackSeq[track_idx].extra != -1)
				continue;

			//Check if we've seen this track 
			int seen = tracks_seen[track_idx];

			if (seen == -1) 
			{
				//We haven't yet seen this track, create a new track
				tracks_seen[track_idx] = (int) new_tracks.size();

				ImageKeyVector track;
				track.push_back(ImageKey(i, j));
				new_tracks.push_back(track);
				track_idxs.push_back(track_idx);

				//printf("%d ", track_idx);
			}
			else 
			{
				new_tracks[seen].push_back(ImageKey(i, j));
			}
		}
	}
	printf("\n");

	delete [] tracks_seen;

	//Now for each (sub) track, triangulate to see if the track is
	
	int pt_count = num_points;

	int num_ill_conditioned = 0;
	int num_high_reprojection = 0;
	int num_cheirality_failed = 0;
	int num_added = 0;

	int num_tracks = (int) new_tracks.size();
	for (int i = 0; i < num_tracks; i++) 
	{
		int num_views = (int) new_tracks[i].size();

		if (num_views < min_views) continue;  /* Not enough views */

		/* Check if at least two cameras fix the position of the point */
		bool conditioned = false;
		bool good_distance = false;
		double max_angle = 0.0;
		for (int j = 0; j < num_views; j++) 
		{
			for (int k = j+1; k < num_views; k++) 
			{
				int camera_idx1 = new_tracks[i][j].first;
				int image_idx1 = added_order[camera_idx1];
				int key_idx1 = new_tracks[i][j].second;

				int camera_idx2 = new_tracks[i][k].first;
				int image_idx2 = added_order[camera_idx2];
				int key_idx2 = new_tracks[i][k].second;

				//Keypoint &key1 = GetKey(image_idx1, key_idx1);
				//Keypoint &key2 = GetKey(image_idx2, key_idx2);
				FeatPoint &key1 = imageData[image_idx1]->GetKeyPoint(key_idx1);//GetKey(image_idx1, key_idx1);
				FeatPoint &key2 = imageData[image_idx2]->GetKeyPoint(key_idx2);//GetKey(image_idx2, key_idx2);

				v2_t p = v2_new(key1.cx, key1.cy);
				v2_t q = v2_new(key2.cx, key2.cy);

				/*
				if (m_optimize_for_fisheye) 
				{
					double p_x = Vx(p), p_y = Vy(p);
					double q_x = Vx(q), q_y = Vy(q);

					m_image_data[image_idx1].
						UndistortPoint(p_x, p_y, Vx(p), Vy(p));
					m_image_data[image_idx2].
						UndistortPoint(q_x, q_y, Vx(q), Vy(q));
				}
				*/

				double angle = ComputeRayAngle(p, q, 
					cameras[camera_idx1], 
					cameras[camera_idx2]);

				if (angle > max_angle)
					max_angle = angle;

				/* Check that the angle between the rays is large
				* enough */
				if (RAD2DEG(angle) >= RAY_ANGLE_THRESHOLD) 
				{
					conditioned = true;
				}
				good_distance = true;
			}
		}

		if (!conditioned || !good_distance) 
		{
			num_ill_conditioned++;
			continue;
		}

		double error;
		v3_t pt;

		//
		bool panorama_mode = false;
		if (!panorama_mode) 
		{
			pt = TriangulateNViews(imageData, new_tracks[i], added_order, cameras, 
				error, true);
		}
		else 
		{
			pt = GeneratePointAtInfinity(imageData, new_tracks[i], added_order, cameras, 
				error, true);
		}

		if (isnan(error) || error > max_reprojection_error) 
		{
			num_high_reprojection++;
			continue;	    
		}

		bool all_in_front = true;
		for (int j = 0; j < num_views; j++) 
		{
			int camera_idx = new_tracks[i][j].first;
			bool in_front = CheckCheirality(pt, cameras[camera_idx]);
			if (!in_front) 
			{
				all_in_front = false;
				break;
			}
		}

		if (!all_in_front) 
		{
			num_cheirality_failed++;

			continue;
		}

	

		fflush(stdout);

		points[pt_count] = pt;

		int camera_idx = new_tracks[i][0].first;
		int image_idx = added_order[camera_idx];
		int key_idx = new_tracks[i][0].second;


		unsigned char r = 255; //GetKey(image_idx, key_idx).m_r;
		unsigned char g = 0;   //GetKey(image_idx, key_idx).m_g;
		unsigned char b = 0;   //GetKey(image_idx, key_idx).m_b;
		colors[pt_count] = v3_new((double) r, (double) g, (double) b);

		pt_views.push_back(new_tracks[i]);

		/* Set the point index on the keys */
		for (int j = 0; j < num_views; j++) 
		{
			int camera_idx = new_tracks[i][j].first;
			int image_idx = added_order[camera_idx];
			int key_idx = new_tracks[i][j].second;			
			imageData[image_idx]->SetFeatExtra(key_idx, pt_count);
			//GetKey(image_idx, key_idx).m_extra = pt_count;
		}

		int track_idx = track_idxs[i];
		//m_track_data[track_idx].m_extra = pt_count;
		trackSeq[track_idx].extra = pt_count;

		assert(pt_count<trackSeq.size());

		pt_count++;
		num_added++;
	}

	printf("[AddAllNewPoints] Added %d new points\n", num_added);
	printf("[AddAllNewPoints] Ill-conditioned tracks: %d\n", 
		num_ill_conditioned);
	printf("[AddAllNewPoints] Bad reprojections: %d\n", num_high_reprojection);
	printf("[AddAllNewPoints] Failed cheirality checks: %d\n", 
		num_cheirality_failed);

	return pt_count;
}

//Find the camera with the most matches to existing points
int FindCameraWithMostMatches(int num_cameras, int *added_order,
							  int &parent_idx, int &max_matches,
							  vector<CImageDataBase*> imageData,
							  const vector<TrackInfo>& trackSeq,
							  const vector<ImageKeyVector> &pt_views) 
{
	max_matches = 0;
	int i_best = -1;
	double top_score = 0.0;
	parent_idx = -1;

	int num_images = imageData.size();//GetNumImages();
	for (int i = 0; i < num_images; i++) 
	{
		//if (m_image_data[i].m_ignore_in_bundle)
		//	continue;
		//if (m_only_bundle_init_focal && !m_image_data[i].m_has_init_focal)
		//	continue;

		/* Check if we added this image already */
		bool added = false;
		for (int j = 0; j < num_cameras; j++) 
		{
			if (added_order[j] == i) 
			{
				added = true;
				break;
			}
		}

		if (added)
			continue;

		int num_existing_matches = 0;
		int parent_idx_best = -1;

		//Find the tracks seen by this image 
		const std::vector<int> &tracks = imageData[i]->GetTrackSeq();//m_image_data[i].m_visible_points;
		int num_tracks = (int) tracks.size();

		for (int j = 0; j < num_tracks; j++) 
		{
			int tr = tracks[j];
			if (trackSeq[tr].extra < 0)
				continue;

			/* This tracks corresponds to a point */
			int pt = trackSeq[tr].extra;
			if ((int) pt_views[pt].size() == 0)
				continue;

			num_existing_matches++;
		}

		if (num_existing_matches > 0)
			printf("  existing_matches[%d] = %d\n", i, num_existing_matches);

		double score = num_existing_matches; 

		if (score > 0.0)
			printf("  score[%d]   = %0.3f\n", i, score);

		if (score > top_score) 
		{
			i_best = i;
			parent_idx = parent_idx_best;
			max_matches = num_existing_matches;
			top_score = score;
		}
	}

	if (parent_idx == -1) 
	{
		printf("Error: parent not found\n");
	}
	return i_best;
}

/* Find all cameras with at least N matches to existing points */
vector<ImagePair> FindCamerasWithNMatches(int match_threshold, 
										  int num_cameras, 
										  int *added_order,
										  vector<CImageDataBase*> imageData,
										  const vector<TrackInfo>& trackSeq,
										  const vector<ImageKeyVector> &pt_views) 
{
	std::vector<ImagePair> image_pairs;

	int num_images = imageData.size(); //GetNumImages();
	for (int i = 0; i < num_images; i++) 
	{
		//if (m_image_data[i].m_ignore_in_bundle)
		//	continue;
		//if (m_only_bundle_init_focal && !m_image_data[i].m_has_init_focal)
		//	continue;

		/* Check if we added this image already */
		bool added = false;
		for (int j = 0; j < num_cameras; j++) 
		{
		
			if (added_order[j] == i) 
			{
				added = true;
				break;
			}
		}

		if (added)
			continue;

		int num_existing_matches = 0;
		int parent_idx_best = -1;

		/* Find the tracks seen by this image */
		const std::vector<int> &tracks = imageData[i]->GetTrackSeq();//m_image_data[i].m_visible_points;
		int num_tracks = (int) tracks.size();

		for (int j = 0; j < num_tracks; j++) 
		{
			int tr = tracks[j];
			if (trackSeq[tr].extra < 0)
				continue;
			/* This tracks corresponds to a point */
			int pt = trackSeq[tr].extra;
			if ((int) pt_views[pt].size() == 0)
				continue;
			num_existing_matches++;
		}

		if (num_existing_matches >= match_threshold)
			image_pairs.push_back(ImagePair(i, parent_idx_best));
	}

	return image_pairs;
}


vector<int> RefineCameraParameters( int num_points, 
								v3_t *points, v2_t *projs, 
								int *pt_idxs, camera_params_t *camera,
								double *error_out, 
								bool adjust_focal,
								bool remove_outliers,
								bool optimize_for_fisheye,
								bool estimate_distortion,
								double min_proj_error_threshold,
								double max_proj_error_threshold)
{
	int num_points_curr = num_points;
	v3_t *points_curr = new v3_t[num_points];
	v2_t *projs_curr = new v2_t[num_points];

	memcpy(points_curr, points, num_points * sizeof(v3_t));
	memcpy(projs_curr, projs, num_points * sizeof(v2_t));

	std::vector<int> inliers;

	for (int i = 0; i < num_points; i++)
		inliers.push_back(i);

	int round = 0;

	// First refine with the focal length fixed 
	camera_refine(num_points_curr, points_curr, projs_curr, camera, 0, 0);    

	while (1)
	{
		//printf("[RefineCameraParameters] Calling with %d points\n", num_points_curr);

		camera_refine(num_points_curr, points_curr, projs_curr, camera, 
			adjust_focal ? 1 : 0, estimate_distortion ? 1 : 0);

		if (!remove_outliers)
			break;

		v3_t *points_next = new v3_t[num_points];
		v2_t *projs_next = new v2_t[num_points];

		int count = 0;
		double error = 0.0;
		std::vector<int> inliers_next;

		double *errors = new double[num_points_curr];

		for (int i = 0; i < num_points_curr; i++) 
		{
			v2_t pr = sfm_project_final(camera, points_curr[i], 1, estimate_distortion ? 1 : 0);

			/*
			if (optimize_for_fisheye)
			{
				// Distort pr
				double x = Vx(pr);
				double y = Vy(pr);
				data.DistortPoint(x, y, Vx(pr), Vy(pr));
			}*/

			double dx = Vx(pr) - Vx(projs_curr[i]);
			double dy = Vy(pr) - Vy(projs_curr[i]);
			double diff = sqrt(dx * dx + dy * dy);

			errors[i] = diff;
			error += diff;
		}

		//printf("[RefineCameraParameters] Error: %0.3f\n", error / num_points_curr);

		// Sort and histogram errors
		double med = kth_element_copy(num_points_curr, 
			iround(0.95 * num_points_curr),
			errors);

		// We will tolerate any match with projection error < 8.0
		double threshold = 1.2 * NUM_STDDEV * med; // k * stddev
		threshold = CLAMP(threshold, min_proj_error_threshold, 
			max_proj_error_threshold);  

		// double threshold = min_proj_error_threshold;
		// double threshold = MAX(8.0, med);
		//printf("[RefineCameraParameters] Threshold = %0.3f\n", threshold);
		for (int i = 0; i < num_points_curr; i++) 
		{
			if (errors[i] < threshold) 
			{
				inliers_next.push_back(inliers[i]);

				points_next[count] = points_curr[i];
				projs_next[count] = projs_curr[i];

				count++;
			}
			/*else 
			{
				if (pt_idxs != NULL)
				{
					printf("[RefineCameraParameters] Removing point [%d] with "
						"reprojection error %0.3f\n", pt_idxs[i], errors[i]);
				}
				else 
				{
					printf("[RefineCameraParameters] Removing point with "
						"reprojection error %0.3f\n", errors[i]);
				}
			}*/
		}

#if 1
		qsort(errors, num_points_curr, sizeof(double), compare_doubles);

		double pr_min = errors[0];
		double pr_max = errors[num_points_curr-1];
		double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

		// Break histogram into 10 bins 
		int idx_count = 0;
		for (int i = 0; i < NUM_ERROR_BINS; i++) {
			double max = pr_min + (i+1) * pr_step;
			int start = idx_count;

			while (idx_count < num_points_curr && errors[idx_count] <= max)
				idx_count++;

			int bin_size = idx_count - start;
			
			//printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n", 
			//	max - pr_step, max, bin_size, 
			//	bin_size / (double) num_points_curr);
		}
#endif

		delete [] points_curr;
		delete [] projs_curr;
		delete [] errors;

		points_curr = points_next;
		projs_curr = projs_next;

		if (count == num_points_curr)
			break;  // We're done 

		num_points_curr = count;

		inliers = inliers_next;

		if (count == 0) // Out of measurements
			break;

		round++;

		if (error_out != NULL) 
		{
			*error_out = error;
		}
	}

	//printf("[RefineCameraParameters] Exiting after %d rounds "
	//	"with %d / %d points\n", round + 1, num_points_curr, num_points);

	delete [] points_curr;
	delete [] projs_curr;

	return inliers;
}


double RefinePoints(int num_points, v3_t *points, v2_t *projs,
					int *pt_idxs, camera_params_t *cameras,
					int *added_order,
					vector<CImageDataBase*> imageData,
					const std::vector<ImageKeyVector> &pt_views,
					camera_params_t *camera_out)
{
	bool estimate_distortion = true;
	double error = 0.0;

	/* Triangulate each of the points */
	for (int i = 0; i < num_points; i++) 
	{
		int pt_idx = pt_idxs[i];

		int num_views = (int) pt_views[pt_idx].size() + 1;

		if (num_views < 2) continue;

		v2_t *pv = new v2_t[num_views];
		double *Rs = new double[9 * num_views];
		double *ts = new double[3 * num_views];

		for (int j = 0; j < num_views; j++) {
			camera_params_t *cam = NULL;

			if (j < num_views - 1) 
			{
				int camera_idx = pt_views[pt_idx][j].first;
				int image_idx = added_order[camera_idx];
				int key_idx = pt_views[pt_idx][j].second;

				 
				FeatPoint fpt = imageData[image_idx]->GetKeyPoint(key_idx); //GetKey(image_idx, key_idx);

				double p3[3] = { fpt.cx, fpt.cy, 1.0 };
				double K[9], Kinv[9];
				GetIntrinsics(cameras[camera_idx], K);
				matrix_invert(3, K, Kinv);

				double p_n[3];
				matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

				pv[j] = v2_new(p_n[0], p_n[1]);

				cam = cameras + camera_idx;
			}
			else 
			{
				double p3[3] = { Vx(projs[i]), Vy(projs[i]), 1.0 };
				double K[9], Kinv[9];
				GetIntrinsics(*camera_out, K);
				matrix_invert(3, K, Kinv);

				double p_n[3];
				matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

				pv[j] = v2_new(p_n[0], p_n[1]);
				cam = camera_out;
			}

			memcpy(Rs + 9 * j, cam->R, 9 * sizeof(double));

			matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * j);
			matrix_scale(3, 1, ts + 3 * j, -1.0, ts + 3 * j);
		}

		// points[i] = triangulate_n(num_views, pv, Rs, ts, &error_curr);
		double error_curr = 0.0;
		points[i] = triangulate_n_refine(points[i], num_views, pv, Rs, ts, 
			&error_curr);

		v2_t pr = sfm_project_final(camera_out, points[i], 1,
			estimate_distortion ? 1 : 0);

		double dx = Vx(pr) - Vx(projs[i]);
		double dy = Vy(pr) - Vy(projs[i]);

		error += dx * dx + dy * dy;

		delete [] pv;
		delete [] Rs;
		delete [] ts;
	}

	return sqrt(error / num_points);
}


vector<int> RefineCameraAndPoints(int num_points,
								  v3_t *points, v2_t *projs,
								  int *pt_idxs, 
								  camera_params_t *cameras,
								  int *added_order,
								  vector<CImageDataBase*> imageData,
								  std::vector<ImageKeyVector> &pt_views,
								  camera_params_t *camera_out,
								  bool remove_outliers)
{
	// double error_thresh = 1.0e-6;
	double error_old = DBL_MAX;
	double derror;

	bool fixed_focal_length = false;
	bool optimize_for_fisheye = false;
	bool estimate_distortion = true;

	bool removed;
	std::vector<int> inliers_out;

	for (int i = 0; i < num_points; i++)
		inliers_out.push_back(i);

	int num_points_curr = num_points;
	v3_t *points_curr = new v3_t[num_points];
	v2_t *projs_curr = new v2_t[num_points];
	int *pt_idxs_curr = new int[num_points];

	memcpy(points_curr, points, sizeof(v3_t) * num_points);
	memcpy(projs_curr, projs, sizeof(v2_t) * num_points);
	memcpy(pt_idxs_curr, pt_idxs, sizeof(int) * num_points);

	do {
		removed = false;
		// do {
		{
			double error = 0.0;

			//Refine the camera
			RefineCameraParameters(num_points_curr, 
				points_curr, projs_curr, 
				pt_idxs_curr, camera_out, &error, 
				!fixed_focal_length, false, // true,
				optimize_for_fisheye,
				estimate_distortion,
				MIN_PROJ_ERROR_THRESHOLD,
				MAX_PROJ_ERROR_THRESHOLD);


#if 1 // Change to 0 to skip point polishing
			//Refine the points
			error = RefinePoints(num_points_curr, points_curr, projs_curr, 
				pt_idxs_curr, cameras, added_order, 
				imageData, pt_views, camera_out);

			printf("[RefineCameraAndPoints] "
				"Error (after point polishing): %0.3f\n", error);

			derror = error_old - error;
			error_old = error;
#endif
		}
		// } while (derror > error_thresh);

		if (remove_outliers)
		{
			//Refine the camera once more and remove outliers
			std::vector<int> inliers;
			inliers = RefineCameraParameters(num_points_curr, 
				points_curr, projs_curr,
				pt_idxs_curr, camera_out, 
				NULL, !fixed_focal_length, true,
				optimize_for_fisheye,
				estimate_distortion,
				MIN_PROJ_ERROR_THRESHOLD,
				MAX_PROJ_ERROR_THRESHOLD);

			int num_inliers = inliers.size();

			std::vector<int> inliers_out_next;
			if (num_inliers < num_points_curr) 
			{
				removed = true;

				for (int i = 0; i < num_inliers; i++) 
				{
					points_curr[i] = points_curr[inliers[i]];
					projs_curr[i] = projs_curr[inliers[i]];
					pt_idxs_curr[i] = pt_idxs_curr[inliers[i]];
					inliers_out_next.push_back(inliers_out[i]);
				}

				num_points_curr = num_inliers;
				inliers_out = inliers_out_next;
			}
		}
		else 
		{
			RefineCameraParameters( num_points_curr, 
				points_curr, projs_curr,
				pt_idxs_curr, camera_out, 
				NULL, !fixed_focal_length, false,
				optimize_for_fisheye,
				estimate_distortion,
				MIN_PROJ_ERROR_THRESHOLD,
				MAX_PROJ_ERROR_THRESHOLD);
		}

	} while (removed);

	delete [] points_curr;
	delete [] projs_curr;
	delete [] pt_idxs_curr;

	return inliers_out;
}


// single image pose estimation using the projections between 3D and 2d 
bool FindAndVerifyCamera(int num_points, v3_t *points_solve, v2_t *projs_solve,
						 int *idxs_solve,
						 double *K, double *R, double *t, 
						 double proj_estimation_threshold,
						 double proj_estimation_threshold_weak,
						 std::vector<int> &inliers,
						 std::vector<int> &inliers_weak,
						 std::vector<int> &outliers)
{
	/* First, find the projection matrix */
	double P[12];
	int r = -1;

	if (num_points >= 9) 
	{
		r = find_projection_3x4_ransac(num_points, 
			points_solve, projs_solve, 
			P, /* 2048 */ 4096 /* 100000 */, 
			proj_estimation_threshold);
	}

	if (r == -1) {
		printf("[FindAndVerifyCamera] Couldn't find projection matrix\n");
		return false;
	}

	/* If number of inliers is too low, fail */
	if (r <= MIN_INLIERS_EST_PROJECTION) {
		printf("[FindAndVerifyCamera] Too few inliers to use "
			"projection matrix\n");
		return false;
	}

	double KRinit[9], Kinit[9], Rinit[9], tinit[3];
	memcpy(KRinit + 0, P + 0, 3 * sizeof(double));
	memcpy(KRinit + 3, P + 4, 3 * sizeof(double));
	memcpy(KRinit + 6, P + 8, 3 * sizeof(double));

	dgerqf_driver(3, 3, KRinit, Kinit, Rinit);	    

	/* We want our intrinsics to have a certain form */
	FixIntrinsics(P, Kinit, Rinit, tinit);
	matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

	printf("[FindAndVerifyCamera] Estimated intrinsics:\n");
	matrix_print(3, 3, Kinit);
	printf("[FindAndVerifyCamera] Estimated extrinsics:\n");
	matrix_print(3, 3, Rinit);
	matrix_print(1, 3, tinit);
	fflush(stdout);

	/* Check cheirality constraint */
	printf("[FindAndVerifyCamera] Checking consistency...\n");

	double Rigid[12] = 
	{ Rinit[0], Rinit[1], Rinit[2], tinit[0],
	Rinit[3], Rinit[4], Rinit[5], tinit[1],
	Rinit[6], Rinit[7], Rinit[8], tinit[2] };

	int num_behind = 0;
	for (int j = 0; j < num_points; j++) {
		double p[4] = { Vx(points_solve[j]), 
			Vy(points_solve[j]),
			Vz(points_solve[j]), 1.0 };
		double q[3], q2[3];

		matrix_product(3, 4, 4, 1, Rigid, p, q);
		matrix_product331(Kinit, q, q2);

		double pimg[2] = { -q2[0] / q2[2], -q2[1] / q2[2] };
		double diff = 
			(pimg[0] - Vx(projs_solve[j])) * 
			(pimg[0] - Vx(projs_solve[j])) + 
			(pimg[1] - Vy(projs_solve[j])) * 
			(pimg[1] - Vy(projs_solve[j]));

		diff = sqrt(diff);

		if (diff < proj_estimation_threshold)
			inliers.push_back(j);

		if (diff < proj_estimation_threshold_weak) {
			inliers_weak.push_back(j);
		} else {
			printf("[FindAndVerifyCamera] Removing point [%d] "
				"(reproj. error = %0.3f)\n", idxs_solve[j], diff);
			outliers.push_back(j);
		}

		// EDIT!!!
		if (q[2] > 0.0)
			num_behind++;  /* Cheirality constraint violated */
	}

	if (num_behind >= 0.9 * num_points) {
		printf("[FindAndVerifyCamera] Error: camera is pointing "
			"away from scene\n");
		return false;
	}

	memcpy(K, Kinit, sizeof(double) * 9);
	memcpy(R, Rinit, sizeof(double) * 9);
	memcpy(t, tinit, sizeof(double) * 3);

	// #define COLIN_HACK
#ifdef COLIN_HACK
	matrix_ident(3, R);
	t[0] = t[1] = t[2] = 0.0;
#endif

	return true;
}

void SetFocalConstraint(double initFocus, camera_params_t *params)
{
	params->constrained[6] = true;
	params->constraints[6] = initFocus;
	params->weights[6] = CONSTRAIN_FOCAL_WEIGHT;
}

camera_params_t BundleInitializeImage( const vector<TrackInfo>& trackSeq,
									   int image_idx, int camera_idx,
									   int num_cameras, int num_points,
								       int *added_order,
								       v3_t *points,
									   camera_params_t *parent,
									   camera_params_t *cameras,
									   vector<CImageDataBase*> imageData,
									   vector<ImageKeyVector> &pt_views,
									   bool *success_out,
									   bool refine_cameras_and_points)
{
	//clock_t start = clock();
	if (success_out != NULL)
		*success_out = true;

	CImageDataBase* pNewImageData = imageData[image_idx];

	bool fixed_focal_length = false;
	bool optimize_for_fisheye = false;
	bool estimate_distortion = true;

	/* Load the keys */
	//data.LoadKeys(false, !m_optimize_for_fisheye);
	//SetTracks(image_idx);

	/* **** Connect the new camera to any existing points **** */
	int num_pts_solve = 0;
	int num_keys = (int) pNewImageData->GetFeatNumber();//data.m_keys.size();
	v3_t *points_solve = new v3_t[num_keys];
	v2_t *projs_solve = new v2_t[num_keys];
	v2_t *projs_solve_orig = new v2_t[num_keys];
	int *idxs_solve = new int[num_keys];
	int *keys_solve = new int[num_keys];
	camera_params_t dummy;
	memset(&dummy, 0, sizeof(camera_params_t));

	printf("[BundleInitializeImage] "
		"Connecting existing matches...\n");

	/* Find the tracks seen by this image */
	vector<int> &tracks =  pNewImageData->GetTrackSeq(); //data.m_visible_points;
    vector<int> &ptIndex = pNewImageData->GetPtSeq();
	int num_tracks = (int) tracks.size();

	for (int i = 0; i < num_tracks; i++) 
	{
		int tr = tracks[i];
		if (trackSeq[tr].extra < 0)
			continue;

		/* This tracks corresponds to a point */
		int pt = trackSeq[tr].extra;
		if ((int) pt_views[pt].size() == 0)
			continue;

		int key =  ptIndex[i]; //data.m_visible_keys[i];
		// printf("  Connecting existing point [%d] (cam: %d)\n", 
		//        pt, image_idx);
		/* Add the point to the set we'll use to solve for
		* the camera position */
		points_solve[num_pts_solve] = points[pt];
		/*
		if (m_optimize_for_fisheye) 
		{
			double x = data.m_keys[key].m_x;
			double y = data.m_keys[key].m_y;
			double x_u, y_u;
			data.UndistortPoint(x, y, x_u, y_u);
			projs_solve[num_pts_solve] = v2_new(x_u, y_u);
			projs_solve_orig[num_pts_solve] = v2_new(x, y);
		}
		else
		*/
		{   
			FeatPoint keypt = pNewImageData->GetKeyPoint(key);
			projs_solve[num_pts_solve] = v2_new(keypt.cx, keypt.cy); //v2_new( data.m_keys[key].m_x, data.m_keys[key].m_y);
		}
		idxs_solve[num_pts_solve] = pt;
		keys_solve[num_pts_solve] = key;
		num_pts_solve++;
	}

	if (num_pts_solve < MIN_MAX_MATCHES) 
	{
		printf("[BundleInitializeImage] Couldn't initialize\n");
		if (success_out != NULL)
			*success_out = false;
		delete [] points_solve;
		delete [] projs_solve;
		delete [] projs_solve_orig;
		delete [] idxs_solve;
		delete [] keys_solve;

		//m_image_data[image_idx].UnloadKeys();
		return dummy;
	}

	/* **** Solve for the camera position **** */
	printf("[BundleInitializeImage] Initializing camera...\n");
	fflush(stdout);
	double Kinit[9], Rinit[9], tinit[3];
	std::vector<int> inliers, inliers_weak, outliers;
	bool success = 
		FindAndVerifyCamera(num_pts_solve, points_solve, projs_solve,
		idxs_solve, Kinit, Rinit, tinit, 
		PROJECTION_ESTIMATION_THRESHOLD, 
		16.0 * PROJECTION_ESTIMATION_THRESHOLD, /*4.0*/
		inliers, inliers_weak, outliers);

	if (!success) 
	{
		printf("[BundleInitializeImage] Couldn't initialize\n");
		if (success_out != NULL)
			*success_out = false;
		//camera_params_t dummy;
		delete [] points_solve;
		delete [] projs_solve;
		delete [] projs_solve_orig;
		delete [] idxs_solve;
		delete [] keys_solve;
		//m_image_data[image_idx].UnloadKeys();
		return dummy;
	}

	camera_params_t camera_new;
	InitializeCameraParams(camera_new);
	SetConstraints(&camera_new, true, true);
	
	/* Start with the new camera at same place as the best match */
	if (success) 
	{
		/* Set up the new camera */
		memcpy(camera_new.R, Rinit, 9 * sizeof(double));
		matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, camera_new.t);
		matrix_scale(3, 1, camera_new.t, -1.0, camera_new.t);
		
		/* Set up the new focal length */
		
		if (pNewImageData->IsHasInitFocus()) 
		{
			double ratio;
			double init = pNewImageData->GetInitFocus() ;//data.m_init_focal;
			double obs = 0.5 * (Kinit[0] + Kinit[4]);

			printf("[BundleInitializeImage] " "Camera has initial focal length of %0.3f\n ", init);

			if (init > obs) ratio = init / obs;
			else            ratio = obs / init;

			if (ratio < 1.4 ) 
			{
				camera_new.f = init; //data.m_init_focal;
				SetFocalConstraint(init, &camera_new);
			}
			else 
			{
				printf("[BundleInitializeImage] "
					"Estimated focal length of %0.3f "
					"is too different\n", obs);
				camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
			}
		}
		else 
		{
			camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
		}
		
		//SetCameraConstraints(added_order[num_cameras], &camera_new);				
		//camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
		/*
		if (m_fixed_focal_length) 
		{
			camera_new.f = m_init_focal_length;
		}
		else
		{
			if (use_focal_estimate) 
			{
				if (data.m_has_init_focal) 
				{
					double ratio;
					double init = data.m_init_focal;
					double obs = 0.5 * (Kinit[0] + Kinit[4]);

					printf("[BundleInitializeImage] "
						"Camera has initial focal length of %0.3f\n", init);

					if (init > obs) ratio = init / obs;
					else            ratio = obs / init;

					if (ratio < 1.4 || m_trust_focal_estimate) 
					{
						camera_new.f = data.m_init_focal;
						if (m_constrain_focal)
							SetFocalConstraint(m_image_data[image_idx], &camera_new);
					}
					else 
					{
						printf("[BundleInitializeImage] "
							"Estimated focal length of %0.3f "
							"is too different\n", obs);
						camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
					}
				}
				else 
				{
					camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
				}
				// } else if (data.m_has_init_focal) {
				//    camera_new.f = data.m_init_focal;
			}
			else 
			{
				if (parent != NULL)
					camera_new.f = parent->f;
				else
					camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
				// 1.2 * MAX(data.GetWidth(), data.GetHeight());
			}
		}*/		
	}
	else 
	{
		printf("[BundleInitializeImage] Error!  "
			"Pose estimation failed!\n");
	}

	/* **** Finally, start the bundle adjustment **** */
	printf("[BundleInitializeImage] Adjusting...\n");
	fflush(stdout);

	int num_inliers = (int) inliers_weak.size();
	v3_t *points_final = new v3_t[num_inliers];
	v2_t *projs_final = new v2_t[num_inliers];
	int *idxs_final = new int[num_inliers];
	int *keys_final = new int[num_inliers];
	int num_points_final = num_inliers;
	for (int i = 0; i < num_inliers; i++) 
	{
		points_final[i] = points_solve[inliers_weak[i]];
		if (false)
			projs_final[i] = projs_solve_orig[inliers_weak[i]];
		else
			projs_final[i] = projs_solve[inliers_weak[i]];
		idxs_final[i] = idxs_solve[inliers_weak[i]];
		keys_final[i] = keys_solve[inliers_weak[i]];
	}

	if (refine_cameras_and_points) 
	{
		inliers = 
			RefineCameraAndPoints(num_points_final,
			points_final, projs_final, idxs_final, 
			cameras, added_order, 
			imageData,
			pt_views, &camera_new, 
			true);
	}
	else
	{
		inliers = RefineCameraParameters(num_points_final, 
			points_final, projs_final, 
			idxs_final, &camera_new, 
			NULL, !fixed_focal_length, true,
			optimize_for_fisheye,
			estimate_distortion,
			MIN_PROJ_ERROR_THRESHOLD,
			MAX_PROJ_ERROR_THRESHOLD);
	}

	if ((int) inliers.size() < 8 || camera_new.f < 0.1 * pNewImageData->GetWd()) 
	{
		printf("[BundleInitializeImage] Bad camera\n");
		if (success_out)
			*success_out = false;
		delete [] points_final;
		delete [] projs_final;
		delete [] idxs_final;
		delete [] keys_final;
		delete [] points_solve;
		delete [] projs_solve;
		delete [] projs_solve_orig;
		delete [] idxs_solve;
		delete [] keys_solve;

		//camera_params_t dummy;
		return dummy;
	}

	/* Point the keys to their corresponding points */
	num_inliers = (int) inliers.size();
	for (int i = 0; i < num_inliers; i++) 
	{
		int inlier_idx = inliers[i];
		// printf("[BundleInitializeImage] Connecting point [%d]\n",
		//        idxs_final[inlier_idx]);
		//data.m_keys[ keys_final[inlier_idx]].m_extra = idxs_final[inlier_idx];
		pNewImageData->SetFeatExtra(keys_final[inlier_idx], idxs_final[inlier_idx]);
		pt_views[idxs_final[inlier_idx]].push_back(ImageKey(camera_idx, keys_final[inlier_idx]));
	}
	fflush(stdout);

	delete [] points_final;
	delete [] projs_final;
	delete [] idxs_final;
	delete [] keys_final;

	delete [] points_solve;
	delete [] projs_solve;
	delete [] projs_solve_orig;
	delete [] idxs_solve;
	delete [] keys_solve;

	//clock_t end = clock();
	//printf("[BundleInitializeImage] Initializing took %0.3fs\n",
	//	(double) (end - start) / CLOCKS_PER_SEC);

	//data.ReadKeyColors();
	//data.m_camera.m_adjusted = true;

	return camera_new;
}


int RemoveBadPointsAndCameras(int num_points, int num_cameras, 
							  int *added_order, 
							  camera_params_t *cameras,
							  v3_t *points, v3_t *colors,
							  vector<CImageDataBase*> imageData,
							  vector<ImageKeyVector> &pt_views)
{
	int num_pruned = 0;

	for (int i = 0; i < num_points; i++) 
	{
		double *pos = points[i].p;
		int num_views = (int) pt_views[i].size();

		if (num_views == 0)
			continue;

		double max_angle = 0.0;
		for (int j = 0; j < num_views; j++) 
		{
			int v1 = pt_views[i][j].first;

			double r1[3];
			matrix_diff(3, 1, 3, 1, pos, cameras[v1].t, r1);
			double norm = matrix_norm(3, 1, r1);
			matrix_scale(3, 1, r1, 1.0 / norm, r1);

			for (int k = j+1; k < num_views; k++) 
			{
				int v2 = pt_views[i][k].first;

				double r2[3];
				matrix_diff(3, 1, 3, 1, pos, cameras[v2].t, r2);
				double norm = matrix_norm(3, 1, r2);
				matrix_scale(3, 1, r2, 1.0 / norm, r2);

				double dot;
				matrix_product(1, 3, 3, 1, r1, r2, &dot);

				double angle = 
					acos(CLAMP(dot, -1.0 + 1.0e-8, 1.0 - 1.0e-8));

				if (angle > max_angle) {
					max_angle = angle;
				}
			}
		}

		if (RAD2DEG(max_angle) < 0.5 * RAY_ANGLE_THRESHOLD) 
		{
			printf("[RemoveBadPointsAndCamera] "
				"Removing point %d with angle %0.3f\n",
				i, RAD2DEG(max_angle));

			for (int j = 0; j < num_views; j++) 
			{
				// Set extra flag back to 0
				int v = pt_views[i][j].first;
				int k = pt_views[i][j].second;
				//GetKey(added_order[v], k).m_extra = -1;
				imageData[added_order[v]]->SetFeatExtra(k,-1);
			}

			pt_views[i].clear();
			if (colors != NULL) 
			{
				Vx(colors[i]) = 0x0;
				Vy(colors[i]) = 0x0;
				Vz(colors[i]) = 0xff;
			}
			num_pruned++;
		}
	}
	printf("[RemoveBadPointsAndCameras] Pruned %d points\n", num_pruned);
	return num_pruned;
}

//bundle adjustment initialization
void InitBundleAdjust(vector<CImageDataBase*> imageData, vector<TrackInfo>& trackSeq)
{
	for(int j=0; j<imageData.size(); j++)
	{
		int numKey = imageData[j]->GetFeatNumber();
		for(int i=0; i<numKey; i++)
		{
			imageData[j]->SetFeatExtra(i, -1);
		}
	}
	for(int i=0; i<trackSeq.size(); i++)
	{
		trackSeq[i].extra = -1;
	}
}

double runSFMApi( int num_pts, int num_cameras, int start_camera,
			   bool fix_points, vector<CImageDataBase*> imageData, 
			   camera_params_t *init_camera_params,
			   v3_t *init_pts, int *added_order, v3_t *colors,
			   std::vector<ImageKeyVector> &pt_views, double eps2, 
			   double *S, double *U, double *V, double *W,
			   bool remove_outliers)
{
	int num_outliers   = 0;
	int total_outliers = 0;
	double dist_total  = 0.0;
	int num_dists      = 0;

	int  *remap   = new int [num_pts];
	v3_t *nz_pts  = new v3_t[num_pts];

	do{
		if (num_pts - total_outliers < MIN_POINTS) 
		{
			printf("[RunSFM] Too few points remaining, exiting!\n");
			fflush(stdout);
			dist_total = DBL_MAX;
			break;
		}

		// Set up the vmask and projections
		char *vmask = NULL;
		double *projections = NULL;

		int num_projections = 0;
		for (int i = 0; i < num_pts; i++) 
		{
			num_projections += (int) pt_views[i].size();
		}

		vmask = new char[num_pts * num_cameras];
		projections = new double[2 * num_projections];

		for (int i = 0; i < num_pts * num_cameras; i++)
			vmask[i] = 0;

		int arr_idx = 0;
		int nz_count = 0;
		for (int i = 0; i < num_pts; i++) 
		{
			int num_views = (int) pt_views[i].size();

			if (num_views > 0) 
			{
				for (int j = 0; j < num_views; j++) 
				{
					int c = pt_views[i][j].first;
					int v = added_order[c];
					int k = pt_views[i][j].second;

					vmask[nz_count * num_cameras + c] = 1;

					projections[2 * arr_idx + 0] = (imageData[v]->GetKeyPoint(k)).cx; //GetKey(v,k).m_x;
					projections[2 * arr_idx + 1] = (imageData[v]->GetKeyPoint(k)).cy; //GetKey(v,k).m_y;

					arr_idx++;
				}

				remap[i] = nz_count;
				nz_pts[nz_count] = init_pts[i];
				nz_count++;
			} 
			else 
			{
				remap[i] = -1;
			}
		}

		dist_total = 0.0;
		num_dists = 0;

		bool fixed_focal = false;
		bool estimate_distortion = true;
        bool use_constraints = false;
		bool constrain_focal = true;
		bool use_point_constraints = false;
		//bool point_constraints = false;
		//bool point_constraint_weight = false;
		bool optimize_for_fisheye = false;
        bool fix_points = false;

		//clock_t start = clock();

		run_sfm(nz_count, num_cameras, start_camera, vmask, projections, 
			    fixed_focal?0:1, 0,
			    estimate_distortion?1:0, 1,
				init_camera_params, nz_pts, 
				(use_constraints || constrain_focal)?1:0,
				(use_point_constraints)?1:0,
				NULL, 
				NULL,
				fix_points?1:0, optimize_for_fisheye, eps2, V, S, U, W);

		//clock_t end = clock();
		//printf("[RunSFM] run_sfm took %0.3fs\n",(double) (end - start) / (double) CLOCKS_PER_SEC);

		// Check for outliers 
		//start = clock();
		std::vector<int> outliers;
		std::vector<double> reproj_errors;
		for (int i = 0; i < num_cameras; i++) 
		{
			int nImageId = added_order[i];

			double K[9] = { init_camera_params[i].f, 0.0, 0.0, 
							0.0, init_camera_params[i].f, 0.0,
							0.0, 0.0, 1.0 };
			double dt[3] = { init_camera_params[i].t[0],
							 init_camera_params[i].t[1],
							 init_camera_params[i].t[2] };		    

			// Compute inverse distortion parameters
			bool bIsEstimate_distortion = true;
			if (bIsEstimate_distortion) 
			{
				double *k = init_camera_params[i].k;
				double k_dist[6] = { 0.0, 1.0, 0.0, k[0], 0.0, k[1] };
				double w_2 = 0.5 * imageData[nImageId]->GetWd(); //data.GetWidth();
				double h_2 = 0.5 * imageData[nImageId]->GetHt(); //data.GetHeight();
				double max_radius = 
					sqrt(w_2 * w_2 + h_2 * h_2) / init_camera_params[i].f;

				InvertDistortion(6, 6, 0.0, max_radius, 
					k_dist, init_camera_params[i].k_inv);
			}

			/*
			bool bIsknown_intrinsics = false;
			if (bIsknown_intrinsics) 
			{
				double *k = init_camera_params[i].k_known;
				double k_dist[8] = 
				{ 0.0, 1.0, 0.0, k[0], 0.0, k[1], 0.0, k[4]};
				double w_2 = 0.5 * data.GetWidth();
				double h_2 = 0.5 * data.GetHeight();
				double max_radius = 
					sqrt(w_2 * w_2 + h_2 * h_2) /
					init_camera_params[i].K_known[0];

				InvertDistortion(8, 6, 0.0, max_radius, k_dist, 
					init_camera_params[i].k_inv);
			}
			*/

			int num_keys =  imageData[ added_order[i] ]->GetFeatNumber(); //GetNumKeys(added_order[i]);

			int num_pts_proj = 0;
			for (int j = 0; j < num_keys; j++) 
			{
				//if (GetKey(added_order[i], j).m_extra >= 0) 
				if(  ( imageData[ added_order[i] ]->GetKeyPoint(j) ).extra>=0 )
				{
					num_pts_proj++;
				}
			}

			double *dists = new double[num_pts_proj];
			int pt_count = 0;

			//std::vector<Keypoint>::iterator iter;
			for (int j = 0; j < num_keys; j++) 
			{
			
				/*for (iter = m_image_data[added_order[i]].m_keys.begin();
				iter != m_image_data[added_order[i]].m_keys.end();
				iter++)
				*/			
				//const Keypoint &key = *iter;
				FeatPoint key = imageData[ added_order[i] ]->GetKeyPoint(j);

				if (key.extra >= 0) 
				{
					double b[3], pr[2];
					double dx, dy, dist;
					int pt_idx = key.extra;

					b[0] = Vx(nz_pts[remap[pt_idx]]);
					b[1] = Vy(nz_pts[remap[pt_idx]]);
					b[2] = Vz(nz_pts[remap[pt_idx]]);

					sfm_project_rd(&(init_camera_params[i]), K, 
						init_camera_params[i].k,
						init_camera_params[i].R, dt, b, pr, 
						estimate_distortion, true);

					if (optimize_for_fisheye) 
					{
						//Distort the points
						double x = pr[0], y = pr[1];
						//m_image_data[added_order[i]].DistortPoint(x, y, pr[0], pr[1]);
					}

					dx = pr[0] - key.cx;
					dy = pr[1] - key.cy;

					dist = sqrt(dx * dx + dy * dy);
					dist_total += dist;
					num_dists++;

					dists[pt_count] = dist;

					pt_count++;
				}
			}

			//Estimate the median of the distances 
			double med = kth_element_copy(num_pts_proj, 
				iround(0.8 * num_pts_proj),
				dists);

			median_copy(num_pts_proj, dists);

			double thresh = 1.2 * NUM_STDDEV * med;  // k * stddev
			thresh = CLAMP(thresh, MIN_PROJ_ERROR_THRESHOLD, MAX_PROJ_ERROR_THRESHOLD);  

			//Compute the average reprojection error for this camera
			double sum = 0.0;
			for (int j = 0; j < num_pts_proj; j++) 
			{
				sum += dists[j];
			}

			double avg = sum / num_pts_proj;
			/*
			printf("[RunSFM] Mean error cam %d[%d] [%d pts]: %0.3e "
				"[med: %0.3e, %0.3e]\n",
				i, added_order[i], num_pts_proj, avg, 
				kth_element_copy(num_pts_proj, 
				iround(0.5 * num_pts_proj), dists), 
				thresh);
				*/

			// printf("Outlier threshold is %0.3f\n", thresh);

			pt_count = 0;
			for (int j = 0; j < num_keys; j++) 
			{
				int pt_idx = ( imageData[ added_order[i] ]->GetKeyPoint(j) ).extra; //GetKey(added_order[i],j).m_extra;

				if (pt_idx < 0)
					continue;

				/*
				// Don't remove constrained points
				if (use_point_constraints && Vx(m_point_constraints[pt_idx]) != 0.0) 
				{
						pt_count++;
						continue;
				}
				*/
				
				if (dists[pt_count] > thresh) 
				{
					// Remove this point from consideration
					bool found = false;
					for (int k = 0; k < (int) outliers.size(); k++) 
					{
						if (outliers[k] == pt_idx) 
						{
							found = true;
							break;
						}
					}
					if (!found) 
					{
						outliers.push_back(pt_idx);
						reproj_errors.push_back(dists[pt_count]);
					}
				}
				pt_count++;
			}


#ifdef  OUTPUT_VERBOSE_STATS

			qsort(dists, num_pts_proj, sizeof(double), compare_doubles);

			double pr_min = dists[0];
			double pr_max = dists[num_pts_proj-1];
			double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

			// Break histogram into 10 bins
			int idx_count = 0;
			for (int i = 0; i < NUM_ERROR_BINS; i++) 
			{
				double max = pr_min + (i+1) * pr_step;
				int start = idx_count;

				while (idx_count < num_pts_proj && dists[idx_count] <= max)
					idx_count++;

				int bin_size = idx_count - start;
				//printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n", 
				//	max - pr_step, max, bin_size, 
				//	bin_size / (double) num_pts_proj);
			}
#endif

			delete [] dists;
		}

		// Remove outlying points
		if (remove_outliers) 
		{
			for (int i = 0; i < (int) outliers.size(); i++) 
			{
				int idx = outliers[i];

				printf("[RunSFM] Removing outlier %d " "(reproj error: %0.3f)\n", idx, reproj_errors[i]);

				if (colors != NULL) 
				{
					Vx(colors[idx]) = 0x0;
					Vy(colors[idx]) = 0x0;
					Vz(colors[idx]) = 0xff;
				}

				int num_views = (int) pt_views[idx].size();

				for (int j = 0; j < num_views; j++) 
				{
					int v = pt_views[idx][j].first;
					int k = pt_views[idx][j].second;

					vmask[idx * num_cameras + v] = 0;

					// Sanity check
					int extra = (imageData[added_order[v]]->GetKeyPoint(k)).extra;
					if (  extra != idx)
						printf("Error!  Entry for (%d,%d) "
						"should be %d, but is %d\n",
						added_order[v], k,
						idx, extra);

					//GetKey(added_order[v], k).m_extra = -2;
					imageData[added_order[v]]->SetFeatExtra(k,-2);
				}
				pt_views[idx].clear();
			}

			num_outliers = outliers.size();
			total_outliers += num_outliers;

			//end = clock();
			//printf("[RunSFM] outlier removal took %0.3fs\n",(double) (end - start) / (double) CLOCKS_PER_SEC);

			printf("[RunSFM] Removing %d outliers\n", num_outliers);
		}

		delete [] vmask;
		delete [] projections;

		for (int i = 0; i < num_pts; i++) 
		{
			if (remap[i] != -1) 
			{
				init_pts[i] = nz_pts[remap[i]];
			}
		}
		if (!remove_outliers) 
			break;
	} while (num_outliers > 0);

	delete [] remap;
	delete [] nz_pts;	
	return dist_total / num_dists;
}



/* Return the intersection of two int vectors */
std::vector<int> GetVectorIntersection(const std::vector<int> &v1,
									   const std::vector<int> &v2)
{
#ifndef WIN32
	__gnu_cxx::hash_set<int> seen;
#else
	stdext::hash_set<int> seen;
#endif

	int v1_size = (int) v1.size();
	int v2_size = (int) v2.size();

	std::vector<int> intersection;

	for (int i = 0; i < v1_size; i++)
		seen.insert(v1[i]);

	for (int i = 0; i < v2_size; i++) {
		if (seen.find(v2[i]) != seen.end())
			intersection.push_back(v2[i]);
	}

	seen.clear();

	return intersection;
}

int SetConstraints(camera_params_t* cam, bool bIsEstimateFocus, bool bIsEstimateDistort)
{
	memset(cam->constrained, 0, sizeof(char)*NUM_CAMERA_PARAMS);
	memset(cam->constraints, 0, sizeof(double)*NUM_CAMERA_PARAMS);
	memset(cam->weights, 0, sizeof(double)*NUM_CAMERA_PARAMS);

	if(bIsEstimateFocus)
	{
		cam->constrained[6] = 1;
		cam->constraints[6] = cam->f;
		cam->weights[6] = 0.0001;
	}
	if(bIsEstimateDistort)
	{
		cam->constrained[7] = 1;
		cam->constrained[8] = 1;

		cam->constraints[7] = 0;
		cam->constraints[8] = 0;

		cam->weights[7] = 100;
		cam->weights[8] = 100;
	}
	return 1;
}


vector<MatchPairIndex> GetMatchList(int imgId1, int imgId2, vector<CImageDataBase*> imageData)
{
	vector<MatchPairIndex> mp;
	
	//get the track sequences corresponding to the image feature points 
	vector<int> &tracks1 = imageData[imgId1]->GetTrackSeq();
	vector<int> &tracks2 = imageData[imgId2]->GetTrackSeq();
	
	//get the intersection of tracks 
	vector<int> isect = GetVectorIntersection(tracks1, tracks2);
	int num_isect = (int) isect.size();
	if (num_isect == 0)
		return mp;

	mp.clear();
	//mp.resize(num_isect);
	printf("GetMatchList.... \n");
	for (int i = 0; i < num_isect; i++) 
	{
		int tr = isect[i];

		vector<int> &pt1 = imageData[imgId1]->GetTrackSeq();
		vector<int>::iterator result = find( pt1.begin( ), pt1.end( ), tr );
		int offset = result - pt1.begin();		
		int k1 = (imageData[imgId1]->GetPtSeq())[offset];

		vector<int> &pt2 = imageData[imgId2]->GetTrackSeq();
		result = find( pt2.begin( ), pt2.end( ), tr );
		offset = result - pt2.begin();
		int k2 = (imageData[imgId2]->GetPtSeq())[offset];


		/*
		//search the position of track
		std::pair<std::vector<int>::const_iterator, 
			std::vector<int>::const_iterator> p;
		const std::vector<int> &pt1 = imageData[imgId1]->GetTrackSeq();
		p = equal_range(pt1.begin(), pt1.end(), tr);
		assert(p.first != p.second);
		int offset = p.first - pt1.begin();
		int k1 = (imageData[imgId1]->GetPtSeq())[offset];

		const std::vector<int> &pt2 = imageData[imgId2]->GetTrackSeq();
		p = equal_range(pt2.begin(), pt2.end(), tr);
		assert(p.first != p.second);
		offset = p.first - pt2.begin();
		int k2 = (imageData[imgId2]->GetPtSeq())[offset];
		*/

		MatchPairIndex match;
		match.l = k1;
		match.r = k2;
		mp.push_back(match);
		printf(".");
	}
    printf("\n");
	return mp;
}


int CeresBA( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, 
			 vector<int> cameraIDOrder,	vector<CameraPara> &cameras)
{

	int  num_pts = trackSeq.size();
	int  num_cameras = cameraIDOrder.size();

	//collect camera parameters
	double* pInteriorParams = new double[3];          //focal length, k1, k2
	pInteriorParams[0] = cameras[0].focus;
	pInteriorParams[1] = 0;
	pInteriorParams[2] = 0;	
	double* pOuterParams = new double[num_cameras*6]; //omiga, phi, kapa, t0,t1,t2
	for(int i=0; i<num_cameras; i++)
	{
		pOuterParams[i*6]   = cameras[i].ax;
		pOuterParams[i*6+1] = cameras[i].ay;
		pOuterParams[i*6+2] = cameras[i].az;
		pOuterParams[i*6+3] = cameras[i].t[0];
		pOuterParams[i*6+4] = cameras[i].t[1];
		pOuterParams[i*6+5] = cameras[i].t[2];
	}

	//collect ground points
	double* grdPt = new double[num_pts*3]; //(double*)malloc( _pts*3*sizeof(double) );
	for(int i=0; i<num_pts; i++)
	{
		grdPt[i*3]   = trackSeq[i].grd.p[0];
		grdPt[i*3+1] = trackSeq[i].grd.p[1];
		grdPt[i*3+2] = trackSeq[i].grd.p[2];
	}
	
	//generate observation points
	int nProjection = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview = trackSeq[i].views.size();
		nProjection += nview;
	}

	vector<int> vecCamIndex;    // camera index for each projection point
	vector<int> vecTrackIndex;  // track point index for each projection point
	vecCamIndex.resize(nProjection);
	vecTrackIndex.resize(nProjection);

	double* projections = new double[nProjection*2]; //(double*)malloc(nProjection*2*sizeof(double));
	int ip = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{   
		int trackIndex = trackSeq[i].extra;
		int nview = trackSeq[i].views.size();
		for(int j=0; j<nview; j++)
		{
			int cameraID = trackSeq[i].views[j].first;
			//int cameraID = cameraIDOrder[index];		
			int ptIndex  = trackSeq[i].views[j].second;

			projections[ip*2]   = imageFeatures[cameraID].featPts[ptIndex].cx; 
			projections[ip*2+1] = imageFeatures[cameraID].featPts[ptIndex].cy; 			
			ip++;

			vecTrackIndex[ip]  = trackIndex;
			vecCamIndex[ip]    = cameraID;
		}
	}

	
	//invoke ceres functions, each time adding one projection





	printf("3D Points: %d   projection number: %d \n", num_pts, nProjection);


	return 0;
}

int DoSFM( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, vector<int> cameraIDOrder,
	vector<CameraPara> &cameras)
{
	//vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews,
	int  num_pts = trackSeq.size();
	int  num_cameras = cameraIDOrder.size();
	v3_t *init_pts = (v3_t*)malloc(sizeof(v3_t)*num_pts);

	for(int i=0; i<num_pts; i++)
	{
		init_pts[i].p[0] = trackSeq[i].grd.p[0];
		init_pts[i].p[1] = trackSeq[i].grd.p[1];
		init_pts[i].p[2] = trackSeq[i].grd.p[2];
	}

	int nProjection = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview = trackSeq[i].views.size();
		nProjection += nview;
	}

	printf("3D Points: %d   projection number: %d \n", num_pts, nProjection);

	//generate mask to indicate whose cameras the 3D point lies? 
	char* vmask = (char*)malloc(num_pts*num_cameras);
	memset(vmask, 0, sizeof(char));
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview = trackSeq[i].views.size();
		for(int j=0; j<nview; j++)
		{
			int index = trackSeq[i].views[j].first;
			vmask[i*num_cameras + index] = 1;
		}
	}

	//generate projection points
	double* projections = (double*)malloc(nProjection*2*sizeof(double));
	int ip = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview = trackSeq[i].views.size();
		for(int j=0; j<nview; j++)
		{
			int index = trackSeq[i].views[j].first;
			int cameraID = cameraIDOrder[index];		
			int ptIndex  = trackSeq[i].views[j].second;

			projections[ip*2]   = imageFeatures[cameraID].featPts[ptIndex].cx; 
			projections[ip*2+1] = imageFeatures[cameraID].featPts[ptIndex].cy; 			
			ip++;
		}
	}

	//invoke sfm-driver
	camera_params_t *pCameras = (camera_params_t *)malloc(sizeof(camera_params_t)*num_cameras);//new camera_params_t[num_cameras];
	for(int i=0; i<num_cameras; i++)
	{
		InitializeCameraParams(pCameras[i]);
		int camId = cameraIDOrder[i];
		pCameras[i].f = cameras[camId].focus;
		pCameras[i].k[0] = cameras[camId].k1;
		pCameras[i].k[1] = cameras[camId].k2;
		memcpy(pCameras[i].R, cameras[camId].R, sizeof(double)*9);
		memcpy(pCameras[i].t, cameras[camId].t, sizeof(double)*3);
	}    

	int    ncons = 0;
	int    est_focal_length = 1; //not fixed
	int    constant_focus = 0;   //focus is not constant  
	int    undistort = 1;        //distortion is considered
	int    explicit_camera_centers = 1;  //R(X-T)
	int    fix_points = 0;               //?
	int    use_constraints = 1;          //?
	int    use_point_constraints = 0;    //?
	double pt_constraint_weight = 0;   
	int    optimize_for_fisheye = 0;
	double eps2 = 1.0e-12;

	//setup constraints
	for(int i=0; i<num_cameras; i++)
	{
		SetConstraints(pCameras+i, est_focal_length, undistort);
	}

	run_sfm(num_pts, 
		num_cameras, 
		ncons, 
		vmask, 
		projections, 
		est_focal_length, 
		constant_focus, 
		undistort, 
		explicit_camera_centers, 
		pCameras, init_pts, 
		use_constraints, 
		use_point_constraints, 
		NULL, 
		pt_constraint_weight, 
		fix_points, 
		optimize_for_fisheye, 
		eps2, NULL, NULL, NULL, NULL
		);

	//output
	for(int i=0; i<num_cameras; i++)
	{
		printf("Camera: %d \n", i);
		printf("Focus: %lf \n", pCameras[i].f);
		printf("translation:  %lf %lf %lf \n", pCameras[i].t[0], pCameras[i].t[1], pCameras[i].t[2]);
		printf("rotation: \n");
		for(int m=0; m<3; m++)
		{
			for(int n=0; n<3; n++)
			{
				printf("%lf ", pCameras[i].R[m*3+n]);
			}
			printf("\n");
		}
		printf("\n\n");
	}

	for(int i=0; i<num_cameras; i++)
	{
		int camId = cameraIDOrder[i];
		cameras[camId].focus = pCameras[i].f;
		cameras[camId].k1 = pCameras[i].k[0];
		cameras[camId].k2 = pCameras[i].k[1];
		memcpy(cameras[camId].R, pCameras[i].R, sizeof(double)*9);
		memcpy(cameras[camId].t, pCameras[i].t, sizeof(double)*3);
	}    
	for(int i=0; i<num_pts; i++)
	{
		trackSeq[i].grd.p[0] = init_pts[i].p[0];
		trackSeq[i].grd.p[1] = init_pts[i].p[1];
		trackSeq[i].grd.p[2] = init_pts[i].p[2];
	}

	free(projections);
	free(vmask);
	free(init_pts);
	free(pCameras);

	return 1;
}

int TrackSeqImageOrder(vector<TrackInfo>& trackSeq, vector<int> cameraIDOrder)
{

	for(int j=0; j<trackSeq.size(); j++)
	{
		int nview = trackSeq[j].views.size();
		for(int i=0; i<nview; i++)
		{
			int id = trackSeq[j].views[i].first;

			int newId = 0;
			for(int k=0; k<cameraIDOrder.size(); k++)
			{
				if( cameraIDOrder[k]==id )
				{
					newId = k;
					break;
				}
			}
			trackSeq[j].views[i].first = newId;
		}
	}


	return 1;
}


CSBA::CSBA()
{
	//m_estimate_focal_length = 1;
	//m_estimate_distortion = 1;
	//m_constrain_focal = 1;
	//m_use_point_constraints = 0;
}

CSBA::~CSBA()
{

}


/*  x,y,z ------ (img1,id1),(img2,id2),....
     pt3                   ptViews
    
	imageFeatures: projection point 
	cameraIDOrder: the original camera index
    cameras:       intrinsic and exterior parameters for cameras
*/
int CSBA::RunSFM(vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
				 vector<ImgFeature> imageFeatures, vector<int> cameraIDOrder,
				 vector<CameraPara> &cameras)
{
    int  num_pts = pt3.size();
	int  num_cameras = cameras.size();
	v3_t *init_pts = (v3_t*)malloc(sizeof(v3_t)*num_pts);

	for(int i=0; i<num_pts; i++)
	{
		init_pts[i].p[0] = pt3[i].p[0];
		init_pts[i].p[1] = pt3[i].p[1];
		init_pts[i].p[2] = pt3[i].p[2];
	}
	
	int nProjection = 0;
	for(int i=0; i<ptViews.size(); i++)
	{
		int nview = ptViews[i].size();
		nProjection += nview;
	}
	printf("3D Points: %d   projection number: %d \n", num_pts, nProjection);
    
	//generate mask to indicate whose cameras the 3D point lies? 
	char* vmask = (char*)malloc(num_pts*num_cameras);
	memset(vmask, 0, sizeof(char));
	for(int i=0; i<ptViews.size(); i++)
	{
		int nview = ptViews[i].size();
		for(int j=0; j<nview; j++)
		{
			int index = ptViews[i][j].first;
			vmask[i*num_cameras + index] = 1;
		}
	}
    
	//generate projection points
	double* projections = (double*)malloc(nProjection*2*sizeof(double));
	int ip = 0;
	for(int i=0; i<ptViews.size(); i++)
	{
		int nview = ptViews[i].size();
		for(int j=0; j<nview; j++)
		{
			int index = ptViews[i][j].first;
			int cameraID = cameraIDOrder[index];
			int ptIndex  = ptViews[i][j].second;

			projections[ip*2]   = imageFeatures[cameraID].featPts[ptIndex].cx; 
			projections[ip*2+1] = imageFeatures[cameraID].featPts[ptIndex].cy; 			
			ip++;
		}
	}

	//invoke sfm-driver
	camera_params_t *pCameras = (camera_params_t *)malloc(sizeof(camera_params_t)*num_cameras);//new camera_params_t[num_cameras];
    for(int i=0; i<num_cameras; i++)
	{
		InitializeCameraParams(pCameras[i]);
		pCameras[i].f = cameras[i].focus;
		memcpy(pCameras[i].R, cameras[i].R, sizeof(double)*9);
		memcpy(pCameras[i].t, cameras[i].t, sizeof(double)*3);
	}    

	int    ncons = 0;
	int    est_focal_length = 1;         //not fixed
	int    constant_focus = 0;           //focus is not constant  
	int    undistort = 1;                //distortion is considered
	int    explicit_camera_centers = 1;  //R(X-T)
	int    fix_points = 0;               //?
    int    use_constraints = 1;          //?
	int    use_point_constraints = 0;    //?
	double pt_constraint_weight = 0;   
    int    optimize_for_fisheye = 0;
	double eps2 = 1.0e-12;

	//setup constraints
	for(int i=0; i<num_cameras; i++)
	{
		SetConstraints(pCameras+i, est_focal_length, undistort);
	}

	run_sfm(num_pts, 
		    num_cameras, 
			ncons, 
			vmask, 
			projections, 
		    est_focal_length, 
		    constant_focus, 
		    undistort, 
		    explicit_camera_centers, 
			pCameras, init_pts, 
		    use_constraints, 
		    use_point_constraints, 
		    NULL, 
		    pt_constraint_weight, 
		    fix_points, 
		    optimize_for_fisheye, 
		    eps2, NULL, NULL, NULL, NULL
			);

	//output
    for(int i=0; i<num_cameras; i++)
	{
		printf("Camera: %d \n", i);
		printf("Focus: %lf \n", pCameras[i].f);
		printf("translation:  %lf %lf %lf \n", pCameras[i].t[0], pCameras[i].t[1], pCameras[i].t[2]);
		printf("rotation: \n");
		for(int m=0; m<3; m++)
		{
			for(int n=0; n<3; n++)
			{
				printf("%lf ", pCameras[i].R[m*3+n]);
			}
			printf("\n");
		}
		printf("\n\n");
	}


	free(projections);
	free(vmask);
	free(init_pts);
	free(pCameras);

	return 1;
}

/*
//ouput 3D model
CModelFileBase* pOutput = new CPlyModel();
//color setting
vector<Point3DDouble> colors;
for(int i=0; i<gpts.size(); i++)
{
Point3DDouble c;
c.p[0]=255; c.p[1]=0; c.p[2]=0;
colors.push_back(c);
}
//add the camera position into the 3D point set
Point3DDouble pos;
pos.p[0]=cameras[leftImageId].t[0];  pos.p[1] = cameras[leftImageId].t[1];  pos.p[2] = cameras[leftImageId].t[2];
gpts.push_back(pos);
pos.p[0]=cameras[rightImageId].t[0];  pos.p[1] = cameras[rightImageId].t[1];  pos.p[2] = cameras[rightImageId].t[2];
gpts.push_back(pos);
Point3DDouble c;
c.p[0]=0; c.p[1]=255; c.p[2]=0;
colors.push_back(c);
colors.push_back(c);
pOutput->Save("d:\\model.ply", gpts,colors);
*/

int CSBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<ImgFeature>& imageFeatures, 
	vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks)
{



	return 0;
}

int CSBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
				       vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{	

	//1. select the initial pair
	double maxInlier = 0;
	int    index = 0;
	for(int i=0; i<pairMatchs.size(); i++)
	{
		if( maxInlier<pairMatchs[i].inlierRatio )
		{
			maxInlier = pairMatchs[i].inlierRatio;
			index = i;
		}
	}
	//for debug
	//index = 0;
	int leftImageId  = pairMatchs[index].lId;
	int rightImageId = pairMatchs[index].rId;
	vector<MatchPairIndex> mp = GetMatchList(leftImageId, rightImageId, imageData);

	//2. relative pose estimation
	//2.1 pose estimation
	vector<Point2DDouble> lpts,rpts;	
	int nMatch = mp.size();
	for(int i=0; i<nMatch; i++)
	{
		int lpIndex = mp[i].l;
		int rpIndex = mp[i].r;

		FeatPoint p1 = imageData[leftImageId]->GetKeyPoint(lpIndex);
		FeatPoint p2 = imageData[rightImageId]->GetKeyPoint(rpIndex);
		
		Point2DDouble pl,pr;
		pl.p[0] = p1.cx;
		pl.p[1] = p1.cy;
		pr.p[0] = p2.cx;
		pr.p[1] = p2.cy;
		lpts.push_back(pl);
		rpts.push_back(pr);
	}

	CRelativePoseBase* pRP = new CEstimatePose5Point();
	pRP->EstimatePose(lpts, rpts, cameras[leftImageId], cameras[rightImageId] );   
	delete pRP;	

	//2.2 triangulation
	CTriangulateBase* pTri = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	vector<double> errorArray;
	pTri->Triangulate(lpts, rpts, cameras[leftImageId], cameras[rightImageId], gpts, errorArray);
	delete pTri;
	
	//3 bundle adjustment
	
	//3.0 initialization
	InitBundleAdjust(imageData, tracks);

	int num_images = cameras.size();
	int *added_order = new int[num_images];
	camera_params_t *pCamera = new camera_params_t[num_images];
	int  max_pts = (int) tracks.size(); // 1243742; /* HACK! */
	v3_t *points = new v3_t[max_pts];
	v3_t *colors = new v3_t[max_pts];
	vector<ImageKeyVector> pt_views;
	int pt_count;
    
	for(int i=0; i<max_pts; i++)
	{
		colors[i].p[0] = 0;
		colors[i].p[1] = 255;
		colors[i].p[2] = 0;
	}
	InitializeCameraParams(cameras[leftImageId], pCamera[0]);
	InitializeCameraParams(cameras[rightImageId],pCamera[1]);    	
    added_order[0] = leftImageId;
	added_order[1] = rightImageId;
	for(int i=0; i<num_images; i++)
	{
		SetConstraints(pCamera+i, true, true);
	}

	//3.1 construct initial 3D points
	int ptIndex = 0;
	for(int i=0; i<gpts.size(); i++)
	{
		int lpIndex = mp[i].l;
		int rpIndex = mp[i].r;
		int ltrack  = imageData[leftImageId]->GetPointTrackIdx(lpIndex);
		int rtrack  = imageData[rightImageId]->GetPointTrackIdx(rpIndex);		

		assert(ltrack==rtrack);

		if( errorArray[i] > PROJECTION_ESTIMATION_THRESHOLD )
			continue;

		points[ptIndex].p[0] = gpts[i].p[0];
		points[ptIndex].p[1] = gpts[i].p[1];
		points[ptIndex].p[2] = gpts[i].p[2];
		
		ImageKeyVector views;
		views.push_back(ImageKey(0, lpIndex));
		views.push_back(ImageKey(1, rpIndex));
		pt_views.push_back(views);
 
		//generate connection between tracks, feature points and adjusting points
		imageData[leftImageId]->SetFeatExtra(lpIndex, ptIndex);
		imageData[rightImageId]->SetFeatExtra(rpIndex, ptIndex);		
		tracks[ltrack].extra = ptIndex;

		//printf("%d \n", ltrack);

		ptIndex++;
	}
	
	runSFMApi(ptIndex, 2, 0, false, imageData, pCamera, points, added_order, NULL, pt_views);
    
	//3.2 add new images
	int round = 0;
	int curr_num_cameras = 2;
	int curr_num_pts = ptIndex;
	int max_matches = 0;
	while (curr_num_cameras < num_images) 
	{
		int parent_idx;
		int max_cam = FindCameraWithMostMatches(curr_num_cameras, added_order, 
												parent_idx, max_matches, 
												imageData, tracks,
												pt_views);

		printf("[SifterApp::BundleAdjust] max_matches = %d\n", max_matches);

		if (max_matches < MIN_MAX_MATCHES)
			break; /* No more connections */

		/* Find all images with 90% of the matches of the maximum */
		std::vector<ImagePair> image_set;
		
		if (false && max_matches < 48) 
		{
			image_set.push_back(ImagePair(max_cam, parent_idx));
		}
		else
		{
			// int nMatches = MIN(100, iround(0.75 /* 0.9 */ * max_matches));
			int nMatches = iround(0.75 /* 0.9 */ * max_matches);
			image_set = FindCamerasWithNMatches(nMatches,curr_num_cameras, 
										added_order, imageData,
										tracks,	pt_views);
		}

		int num_added_images = (int) image_set.size();
		printf("[SifterApp::BundleAdjustFast] Registering %d images\n", num_added_images);
		for (int i = 0; i < num_added_images; i++)
			printf("[SifterApp::BundleAdjustFast] Adjusting camera %d\n",image_set[i].first);

		/* Now, throw the new cameras into the mix */
		int image_count = 0;
		for (int i = 0; i < num_added_images; i++) 
		{
			int next_idx = image_set[i].first;
			int parent_idx = image_set[i].second;

			added_order[curr_num_cameras + image_count] = next_idx;

			printf("[SifterApp::BundleAdjust[%d]] Adjusting camera %d "
				"(parent = %d)\n", 
				round, next_idx, 
				(parent_idx == -1 ? -1 : added_order[parent_idx]));

			/* **** Set up the new camera **** */
			bool success = false;
			camera_params_t camera_new = 
				BundleInitializeImage(tracks,
				next_idx, curr_num_cameras + image_count,
				curr_num_cameras, curr_num_pts,
				added_order, points, 
				NULL /*cameras + parent_idx*/, 
				pCamera, 
				imageData,
				pt_views, &success);

			if (success) 
			{
				pCamera[curr_num_cameras+image_count] = camera_new;
				//set constraints for new camera, added by Donghai Xie, 2014.1.10
				//SetConstraints(pCamera+curr_num_cameras+image_count, true, true);
				image_count++;
			} 
			else 
			{
				printf("[BundleAdjust] Couldn't initialize image %d\n",	next_idx);
				//m_image_data[next_idx].m_ignore_in_bundle = true;
			}
		}
		
		curr_num_cameras += image_count;

		double dist0 = 0.0;
		pt_count = BundleAdjustAddAllNewPoints( curr_num_pts, curr_num_cameras,
												added_order, pCamera, 
												points, colors,
												dist0, 
												imageData,
												tracks,
												pt_views);	

		curr_num_pts = pt_count;

		/* Run sfm again to update parameters */
		runSFMApi(curr_num_pts, curr_num_cameras, 0, false, imageData,
			pCamera, points, added_order, colors, pt_views);

		/* Remove bad points and cameras */
		RemoveBadPointsAndCameras(curr_num_pts, curr_num_cameras + 1, 
			added_order, pCamera, points, colors, 
			imageData, pt_views);

		//printf("  focal lengths:\n");
		round++;
	}

	//output 	
	char buf[256];
	//sprintf(buf, "points%03d.ply", curr_num_cameras);
	sprintf(buf, "points.ply");
	DumpPointsToPly(outDir, buf, 
		curr_num_pts, curr_num_cameras, 
		points, colors, pCamera);

	//sprintf(buf, "%03d.out", curr_num_cameras);
	sprintf(buf, "bundle.out");
	DumpOutputFile(outDir, buf, 
		num_images, curr_num_cameras, curr_num_pts,
		added_order, pCamera, points, colors, pt_views);

	

	free(points);
	free(colors);
	free(pCamera);
	free(added_order);

	return 1;
}

int CSBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
					   vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{
	//1. select the initial pair
	double maxInlier = 0;
	int    index = 0;
	for(int i=0; i<pairMatchs.size(); i++)
	{
		if( maxInlier<pairMatchs[i].inlierRatio )
		{
			maxInlier = pairMatchs[i].inlierRatio;
			index = i;
		}
	}

	//for debug
	index = 0;

	int leftImageId  = pairMatchs[index].lId;
	int rightImageId = pairMatchs[index].rId;
	vector<TrackInfo> trackSeq;
    GetMatch(leftImageId, rightImageId, tracks, trackSeq);    

	//2. relative pose estimation
	//2.1 pose estimation
	vector<Point2DDouble> lpts,rpts;	
	int nMatch = trackSeq.size();
	for(int i=0; i<nMatch; i++)
	{
		int imgId1 = trackSeq[i].views[0].first;
		int ptId1  = trackSeq[i].views[0].second;
		int imgId2 = trackSeq[i].views[1].first;
		int ptId2  = trackSeq[i].views[1].second;		
		Point2DDouble pl,pr;

		if( imgId1==leftImageId && imgId2==rightImageId )
		{
			pl.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pl.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pr.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			rpts.push_back(pr);
		}
		else
		{
			pl.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pl.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pr.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			rpts.push_back(pr);		
		}
	}

	CRelativePoseBase* pRP = new CEstimatePose5Point();
	pRP->EstimatePose(lpts, rpts, cameras[leftImageId], cameras[rightImageId] );   
		
	//2.2 triangulation
	CTriangulateBase* pTri = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	pTri->Triangulate(lpts, rpts, cameras[leftImageId], cameras[rightImageId], gpts);
    
	//2.3 bundle adjust initialization
	vector<int> cameraIDOrder;
	cameraIDOrder.push_back(leftImageId);
	cameraIDOrder.push_back(rightImageId);
	for(int i=0; i<gpts.size(); i++)
	{
		trackSeq[i].grd = gpts[i];
	}
	DoSFM(trackSeq, imageFeatures, cameraIDOrder, cameras);	
    
	//3. add new cameras and bundle adjustment
	for(int i=0; i<cameras.size()-2; i++)
	{
		//3.1 find the new camera
		int nNewCamId = FindNewImage(cameras.size(), cameraIDOrder, tracks, trackSeq);
		printf("New Camera: %d \n", nNewCamId);

		//3.2 get the 3D and 2D for new image, for DLT 
		vector<Point3DDouble> pt3;
		vector<Point2DDouble> pt2;
		int newId = cameraIDOrder.size();
		UpdateTracks(nNewCamId, newId, tracks, trackSeq, imageFeatures, pt3, pt2);

		//3.3 camera pose estimation
		CPoseEstimationBase* pPoseEstimate = new CDLTPose();
		pPoseEstimate->EstimatePose(pt3, pt2, cameras[nNewCamId]);

		//3.4 add new tracks
		AddNewTracks(nNewCamId, cameraIDOrder, cameras, imageFeatures, tracks, trackSeq);    
		cameraIDOrder.push_back(nNewCamId);

		//3.5 bundle adjustment
		DoSFM(trackSeq, imageFeatures, cameraIDOrder, cameras);		
	}

	return 1;
}


///////////////////////////////////////////////////////////////////////////////////////
CCeresBA::CCeresBA()
{

}

CCeresBA::~CCeresBA()
{

}

int CCeresBA::RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
	vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
	vector<CameraPara>& cameras)
{


	return 0;
}

//
int CCeresBA::BundleAdjust( int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
	vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{


	return 0;
}


//new interface, including more input parameters
int CCeresBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
	vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{


	return 0;
}


int CCeresBA::BundleAdjust(int numCameras, 
							vector<CameraPara>&   cameras,     
							vector<ImgFeature>&   imageFeatures, //feature points
							vector<PairMatchRes>& pairMatchs,    //matching results
							vector<TrackInfo>& tracks)           //track points ( connecting points )
{
	//1. select the initial pair
	double maxInlier = 0;
	int    index = 0;
	for(int i=0; i<pairMatchs.size(); i++)
	{
		if( maxInlier<pairMatchs[i].matchs.size() )
		{
			maxInlier = pairMatchs[i].matchs.size();
			index = i;
		}
	}

	int leftImageId  = pairMatchs[index].lId;
	int rightImageId = pairMatchs[index].rId;
	vector<TrackInfo> trackSeq;
	GetMatch(leftImageId, rightImageId, tracks, trackSeq);    

	//2. relative pose estimation
	//2.1 pose estimation
	vector<Point2DDouble> lpts,rpts;	
	int nMatch = trackSeq.size();
	for(int i=0; i<nMatch; i++)
	{
		int imgId1 = trackSeq[i].views[0].first;
		int ptId1  = trackSeq[i].views[0].second;
		int imgId2 = trackSeq[i].views[1].first;
		int ptId2  = trackSeq[i].views[1].second;		
		Point2DDouble pl,pr;

		if( imgId1==leftImageId && imgId2==rightImageId )
		{
			pl.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pl.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pr.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			rpts.push_back(pr);
		}
		else
		{
			pl.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pl.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pr.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			rpts.push_back(pr);		
		}
	}

	CRelativePoseBase* pRP = new CEstimatePose5Point();
	pRP->EstimatePose(lpts, rpts, cameras[leftImageId], cameras[rightImageId] );   

	//2.2 triangulation
	CTriangulateBase* pTri = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	pTri->Triangulate(lpts, rpts, cameras[leftImageId], cameras[rightImageId], gpts);


	//2.3 bundle adjust for pair
	vector<int> cameraIDOrder;
	cameraIDOrder.push_back(leftImageId);
	cameraIDOrder.push_back(rightImageId);
	for(int i=0; i<gpts.size(); i++)
	{
		trackSeq[i].grd = gpts[i]; //save the initial ground points
	}
	WritePMVSPly("c:\\temp\\pair.ply", gpts);
	CeresBA(trackSeq, imageFeatures, cameraIDOrder, cameras);
	

	return 0;
}