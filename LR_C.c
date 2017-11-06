#include <stdio.h>
#include <stdlib.h>
#include <mat.h>
#include <mkl.h>
#include "lr.h"
#include <sys/time.h>
#include <string.h>

int main()
{
	MATFile *mat_truth = NULL;
	MATFile *mat_D_hyper = NULL;		// MAT file pointer
	MATFile *mat_multi = NULL;
	MATFile *mat_subsample = NULL;
	mxArray *mx_truth = NULL;
	mxArray *mx_D_hyper = NULL;		// mxArray pointer
	mxArray *mx_multi = NULL;
	mxArray *mx_subsample = NULL;

	/******************** CAVE **************************************/
	const char *names[] = { "balloons", "beads", "cd", "chart_and_stuffed_toy",
							"clay", "cloth", "egyptian_statue", "fake_and_real_beers",
							"fake_and_real_food", "fake_and_real_lemon_slices",
							"fake_and_real_lemons", "fake_and_real_peppers",
							"fake_and_real_strawberries", "fake_and_real_sushi",
							"fake_and_real_tomatoes", "feathers", "flowers",
							"glass_tiles", "hairs", "jelly_beans", "oil_painting",
							"paints", "photo_and_face", "pompoms", "real_and_fake_apples",
							"sponges", "stuffed_toys", "superballs", "thread_spools", "watercolors" };

	mat_subsample = matOpen("./original_data_CAVE/subsample.mat", "r");
	/**************************************************************/

	/******************** Havard ***********************************/
	const char *names[] = { "img1", "img2" , "img3", "img4", "img5", "img6",
							"imga1", "imga2", "imga3", "imga4", "imga5", "imga6", "imga7", "imga8",
							"imgb0", "imgb1", "imgb2", "imgb3", "imgb4", "imgb5", "imgb6", "imgb7", "imgb8", "imgb9",
							"imgc1", "imgc2", "imgc3", "imgc4", "imgc5", "imgc6", "imgc7", "imgc8", "imgc9",
							"imgd0", "imgd1", "imgd2", "imgd3", "imgd4", "imgd5", "imgd6", "imgd7", "imgd8", "imgd9",
							"imge0", "imge1", "imge2", "imge3", "imge4", "imge5", "imge6", "imge7",
							"imgf1", "imgf2", "imgf3", "imgf4", "imgf5", "imgf6", "imgf7", "imgf8",
							"imgg0", "imgg1", "imgg2", "imgg3", "imgg4", "imgg5", "imgg6", "imgg7", "imgg8", "imgg9",
							"imgh0", "imgh1", "imgh2", "imgh3", "imgh4", "imgh5", "imgh6", "imgh7" };

	mat_subsample = matOpen("./original_data_Havard/subsample.mat", "r");
	/*************************************************************/

	/****************************** ICVL **************************************/
	const char *names[] = { "4cam_0411-1640-1", "4cam_0411-1648", "BGU_0403-1419-1",
							"BGU_0522-1113-1", "BGU_0522-1127", "BGU_0522-1136", "BGU_0522-1201",
							"BGU_0522-1203", "BGU_0522-1211", "BGU_0522-1216", "BGU_0522-1217",
							"Flower_0325-1336", "Master20150112_f2_colorchecker", "Maz0326-1038",
							"Ramot0325-1364", "bguCAMP_0514-1659", "bguCAMP_0514-1711",
							"bguCAMP_0514-1712", "bguCAMP_0514-1718", "bguCAMP_0514-1723",
							"bguCAMP_0514-1724", "bgu_0403-1439", "bgu_0403-1444", "bgu_0403-1459",
							"bgu_0403-1511", "bgu_0403-1523", "bgu_0403-1525", "eve_0331-1549",
							"eve_0331-1551", "eve_0331-1601", "eve_0331-1602", "eve_0331-1606",
							"eve_0331-1618", "eve_0331-1632", "eve_0331-1633", "eve_0331-1646",
							"eve_0331-1647", "eve_0331-1656", "eve_0331-1657", "eve_0331-1702",
							"eve_0331-1705", "grf_0328-0949", "hill_0325-1219", "hill_0325-1228",
							"hill_0325-1235", "hill_0325-1242", "lst_0408-0950", "lst_0408-1004",
							"lst_0408-1012", "maz_0326-1048", "mor_0328-1209-2", "omer_0331-1055",
							"omer_0331-1102", "omer_0331-1104", "omer_0331-1118", "omer_0331-1119",
							"omer_0331-1130", "omer_0331-1131", "omer_0331-1135", "omer_0331-1150",
							"omer_0331-1159", "pepper_0503-1228", "pepper_0503-1229", "pepper_0503-1236",
							"peppers_0503-1308", "peppers_0503-1311", "peppers_0503-1315",
							"peppers_0503-1330", "peppers_0503-1332", "plt_0411-1037", "plt_0411-1046",
							"plt_0411-1116", "plt_0411-1155", "plt_0411-1200-1", "plt_0411-1207",
							"plt_0411-1210", "plt_0411-1211", "plt_0411-1232-1", "prk_0328-0945",
							"prk_0328-1025", "prk_0328-1031", "prk_0328-1034", "prk_0328-1037",
							"prk_0328-1045", "ramot_0325-1322", "rsh2_0406-1505", "rsh_0406-1343",
							"rsh_0406-1356", "rsh_0406-1413", "rsh_0406-1427", "rsh_0406-1441-1",
							"rsh_0406-1443", "sami_0331-1019", "sat_0406-1107", "sat_0406-1129",
							"sat_0406-1130", "sat_0406-1157-1", "strt_0331-1027" };

	mat_subsample = matOpen("./original_data_ICVL/subsample.mat", "r");
	/*************************************************************************/

	if (mat_subsample == NULL) {
		printf("Error opening MAT file.\n");
		return 1;
	}
	mx_subsample = matGetVariable(mat_subsample, "D");
	matClose(mat_subsample);

	double avr_rmse = 0;
	double avr_time = 0;

	int N = 30;
	for (int i = 0; i < N; i++) {
		/************* CAVE *************/
		char str1[100] = "./original_data_CAVE/";
		char str2[100] = "./original_data_CAVE/";
		char str3[100] = "./original_data_CAVE/";
		/*******************************/

		/*************** Havard ************/
		char str1[100] = "./original_data_Havard/";
		char str2[100] = "./original_data_Havard/";
		char str3[100] = "./original_data_Havard/";
		/*********************************/

		/*************** ICVL ************/
		char str1[100] = "./original_data_ICVL/";
		char str2[100] = "./original_data_ICVL/";
		char str3[100] = "./original_data_ICVL/";
		/*********************************/

		mat_truth = matOpen(strcat(strcat(str1, names[i]), ".mat"), "r");
		mat_D_hyper = matOpen(strcat(strcat(strcat(str2, "D_"), names[i]), ".mat"), "r");
		mat_multi = matOpen(strcat(strcat(strcat(str3, "RGB_"), names[i]), ".mat"), "r");
		if (mat_truth == NULL || mat_D_hyper == NULL || mat_multi == NULL) {
			printf("Error opening MAT file.\n");
			return 1;
		}

		/****************** CAVE **************************/
		mx_truth = matGetVariable(mat_truth, "GroundTruth");
		mx_D_hyper = matGetVariable(mat_D_hyper, "D_GroundTruth");
		mx_multi = matGetVariable(mat_multi, "RGB_GroundTruth");
		/*************************************************/

		/******************** Havard ************************/
		mx_truth = matGetVariable(mat_truth, "ref2");
		mx_D_hyper = matGetVariable(mat_D_hyper, "D_ref2");
		mx_multi = matGetVariable(mat_multi, "RGB_ref2");
		/***************************************************/

		/******************** ICVL ************************/
		mx_truth = matGetVariable(mat_truth, "rad2");
		mx_D_hyper = matGetVariable(mat_D_hyper, "D_rad2");
		mx_multi = matGetVariable(mat_multi, "RGB_rad2");
		/*************************************************/
		
		matClose(mat_truth);
		matClose(mat_multi);
		matClose(mat_D_hyper);

		if (mx_truth == NULL || mx_D_hyper == NULL || mx_multi == NULL || mx_subsample == NULL) {
			printf("mxArray not found!\n");
			return 1;
		}

		const mwSize *size_multi = mxGetDimensions(mx_multi);
		double *estimated_hyper = (double *)mkl_malloc(size_multi[0] * size_multi[1] * 31 * sizeof(double), 32);
		if (estimated_hyper == NULL) {
			printf("Error: can not allocate memory for estimated_hyper! Aborting ......\n");
			return 1;
		}

		struct timeval start, stop;
		gettimeofday(&start, 0);
		hss_lr((const mxArray *)mx_D_hyper, (const mxArray *)mx_multi, (const mxArray *)mx_subsample, estimated_hyper);
		gettimeofday(&stop, 0);

		double RMSE = rmse(mx_truth, estimated_hyper, size_multi[0], size_multi[1], 31);
		avr_rmse += RMSE;
		avr_time += (1000*(stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec)/1000.0);	// ms
		printf("RMSE: %f   run time: %f\n", RMSE, 1000*(stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec)/1000.0);

		mkl_free(estimated_hyper);
		mxDestroyArray(mx_truth);
		mxDestroyArray(mx_multi);
		mxDestroyArray(mx_D_hyper);

		//MATFile *f = matOpen("test.mat", "w");
		//const size_t ndims[3] = { 512, 512, 31 };
		//mxArray *mx_estimated = mxCreateNumericArray(3, ndims, mxDOUBLE_CLASS, mxREAL);
		//double *mx = (double *)mxMalloc(512 * 512 * 31 * sizeof(double));
		//for (int i = 0; i < 512 * 512 * 31; i++) {
		//	mx[i] = estimated_hyper[i];
		//}
		//mxSetPr(mx_estimated, mx);
		//matPutVariable(f, "A", mx_estimated);
		//matClose(f);
	}
	mxDestroyArray(mx_subsample);

	printf("average run time: %f\n", avr_time / N);

    return 0;
}

