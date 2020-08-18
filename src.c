/**************************************************************************
* Implementing k-newMeans clustering on a given dataset of particle positions
*
* @source src.c
* @param
* @author Sharath HP
* @date 10th Nov 2019
* @version 1.6 MPI execution working
* @@
*
**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <dirent.h>
#include <limits.h>
#include <float.h>

#define DEBUG 1

#define ROOT (myRank==0)

#define KMIN 45
#define KMAX 60
#define KSTEP 5

#define D1 17
#define D2 16

/* Points are in 3D */
typedef struct point {
    double x, y, z; // Co-ordinates of the point
	double info; //For point, stands for current cluster index, for mean, stands for # points in the cluster
} point;

/* Compute the squared euclidean distance between given two points */
double distSq ( point p1, point p2 ) {
	double d = 0, delta;
	delta = p1.x - p2.x;
	d += (delta * delta);
	delta = p1.y - p2.y;
	d += (delta * delta);
	delta = p1.z - p2.z;
	d += (delta * delta);
	return d;
}

/* Assign data point to one of k clusters and simultaneously update 
the cluster means partial sum (of co-ordinates) and count (# of process points in this cluster) */
double clusterPoint ( point *data, int k, point *means, point *pMeans ) {
	int index, i;
	double dist, distMin;
	distMin = distSq ( *data, means[0] );
	index = 0;
	for ( i = 1; i < k; i++ ) {
		dist = distSq ( *data, means[i] );
		if ( dist < distMin ) {
			distMin = dist;
			index = i;
		}
	}
	/* Increase no of points in this cluster by 1 */
	pMeans[index].info += 1;
	/* Find cumulative sum of the respective co-ordinates */
	pMeans[index].x += data->x;
	pMeans[index].y += data->y;
	pMeans[index].z += data->z;
	data->info = index;
	return distMin;
}

/* For custom operation in MPI_Reduce */
void combineMeans ( point *pMeans, point *means, int *k, MPI_Datatype *dtype  ) {
	int i;
    for ( i = 0; i < *k; i++ ) {
        means[i].x += pMeans[i].x;
        means[i].y += pMeans[i].y;
        means[i].z += pMeans[i].z;
        means[i].info += pMeans[i].info;
    }
}

/* Zero out the given structure of points */
void zeroMeans ( int n, point *means ) {
	int i;
	for ( i = 0; i < n; i++ ) {
		means[i].info = 0;
		means[i].x = 0;
		means[i].y = 0;
		means[i].z = 0;
	}
}

/* Helper function */
void print ( char *str, int n, point *means ) {
	int i;
	printf ( "\n%s\n", str );
	for ( i = 0; i < n; i++ ) {
		printf ( "Mean %d: %lf %lf %lf\n", (int) means[i].info, means[i].x, means[i].y, means[i].z );
	}
}

void distMin ( point pPoint, point *means, int newK, double *minDist ){
    double dist;
	*minDist = DBL_MAX;
	int i;
    for ( i = 0 ; i < newK; i++ ) {
        dist = distSq (pPoint, means[i]);
        if ( dist < *minDist ) {
            *minDist = dist;
        }
    }
}

void kMeansPP (point *pPoints, int nPoints, point *means, int k ){
    int init = rand() % nPoints;
	/* Initializing first centroid randomly */
    means[0] = pPoints[init];
    int newK = 0, i, j;
	double minDist, dist;

    for ( i = 1; i < k ; i++ ){
        point mean;
        newK++;
        dist = 0;

        for ( j = 0; j < nPoints ; j++ ){
            distMin ( pPoints[j], means, newK, &minDist);
            if (minDist > dist) {
                dist = minDist;
                mean = pPoints[j];
            }
        }
    means[i] = mean;
    }
}

double nTotalTime ( int *kNosRecvCounts, int nProc, double *time, int tStepTotal) {
	int i = 0, j, k;

	double sum = 0;
    double pSum[nProc];
    memset( pSum, 0, sizeof(double) * nProc);

	for ( j = 0; j < nProc; j++ ) {
		for ( k = kNosRecvCounts[j]; k > 0 ; k-- )
			pSum[j] += time[i++];
	}

	double max = 0;
	for ( i = 0; i < nProc; i++ ) {
		if ( pSum[i] > max )
			max = pSum[i];
	}

	return max;
}

double avgTime ( int *kNosRecvCounts, int nProc, double *time, int tStepTotal) {
	int i = 0, j, k;

	double sum = 0, wtSum = 0, wt = 0;
    double pSum[nProc];
    memset( pSum, 0, sizeof(double) * nProc);

	for ( j = 0; j < nProc; j++ ) {
		for ( k = kNosRecvCounts[j]; k > 0 ; k-- )
			pSum[j] += time[i++];
	}

	for ( j = 0; j < nProc; j++ ) {
		wtSum += kNosRecvCounts[j] * pSum[j];
		wt += kNosRecvCounts[j];
	}
	
	wtSum /= wt;

	return wtSum / tStepTotal;
}

int main ( int argc, char *argv[] ) {
	
	/* p - procoess, n - number of, d - data, t - time */
	int myRank, nProc, nPoints, k, error, index, i, iters, n, meanInd;
	int tOffset, tStep, tStepTotal, tpProc;
	int DOUBLE = sizeof (double);
	int POINT = sizeof (point);
	double *pdBuf, threshold;
	double avgPPTime, avgPTime, pTotalTime;
	
	/* Processing points */
	point *pPoints;
	char path[50], kpath[50], run[10], dPath[50], outFiles[3][50], outFolder[50];
	/* Various I/O file handles */
	FILE *fptr, *kfptr, *pptFptr, *ptFptr, *ttFptr;
	/* Optimized values of k, for each of the clusters - pre-computed to preserve same parameters across running different configurations */
	int d1kSteps[] = { 35, 38, 38, 36, 38, 36, 34, 35, 36, 36, 38, 37, 38, 40, 41, 41, 40 };
	int d2kSteps[] = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 };

	/* d - data, f - file */
	MPI_File dFile;
	MPI_Offset dfSize;

    /* MPI initialization, find process rank and total no. of processes */
	MPI_Init ( &argc, &argv );
	MPI_Comm_rank ( MPI_COMM_WORLD, &myRank );
	MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
	
	/* Total Time Start */
	double _totalTime = MPI_Wtime();

	/* Pre Processing Time Start */
	double _timePP = MPI_Wtime();
	
	/* Command Line Argument 1 - total time steps (files) */
	tStepTotal = atoi (argv[1]);
	
	/* Command Line Argument 2 - path to data dump */
	strcpy ( dPath, argv[2] );
	
	/* Command Line Argument 3 - cluster name - CSE/ HPC2010 */
	strcpy ( outFolder, argv[3] );

	/* Output file creation */
	if ( tStepTotal == D1 )
		sprintf ( path, "./%s/data1/output_%d.txt", outFolder, nProc);
	else
		sprintf ( path, "./%s/data2/output_%d.txt", outFolder, nProc);

	if ( tStepTotal == D1 ) {
		sprintf ( outFiles[0], "./%s/data1/ppt.csv", outFolder );
		sprintf ( outFiles[1], "./%s/data1/pt.csv", outFolder );
		sprintf ( outFiles[2], "./%s/data1/tt.csv", outFolder );
	}
	else {
		sprintf ( outFiles[0], "./%s/data2/ppt.csv", outFolder );
		sprintf ( outFiles[1], "./%s/data2/pt.csv", outFolder );
		sprintf ( outFiles[2], "./%s/data2/tt.csv", outFolder );
	}

	if ROOT {
		/* One time write */
		fptr = fopen ( path, "w" );
		if( !fptr ) {
            printf ( "Error Creating %s\n", path );
            exit ( 1 );
        }
		/* Writes across multiple runs */
		if ( access( outFiles[0], F_OK ) != -1 ) {
			pptFptr = fopen ( outFiles[0], "a" );
			if( !pptFptr ) {
				printf ( "Error Opening %s\n", outFiles[0] );
				exit ( 1 );
			}
		}
		else {
			pptFptr = fopen ( outFiles[0], "w" );
			if( !pptFptr ) {
				printf ( "Error Opening %s\n", outFiles[0] );
				exit ( 1 );
			}
			fprintf ( pptFptr, "p,time\n" );
		}
		
		if ( access( outFiles[1], F_OK ) != -1 ) {
			ptFptr = fopen ( outFiles[1], "a" );
			if( !ptFptr ) {
				printf ( "Error Opening %s\n", outFiles[1] );
				exit ( 1 );
			}
		}
		else {
			ptFptr = fopen ( outFiles[1], "w" );
			if( !ptFptr ) {
				printf ( "Error Opening %s\n", outFiles[1] );
				exit ( 1 );
			}
			fprintf ( ptFptr, "p,time\n" );
		}

	
		if ( access( outFiles[2], F_OK ) != -1 ) {
			ttFptr = fopen ( outFiles[2], "a" );
			if( !ttFptr ) {
				printf ( "Error Opening %s\n", outFiles[2] );
				exit ( 1 );
			}
		}
		else {
			ttFptr = fopen ( outFiles[2], "w" );
			if( !ttFptr ) {
				printf ( "Error Opening %s\n", outFiles[2] );
				exit ( 1 );
			}
			fprintf ( ttFptr, "p,time\n" );
		}

	}

	/* Creating a custom data type */
	MPI_Datatype point_t;
	MPI_Type_contiguous ( 4, MPI_DOUBLE, &point_t );
	MPI_Type_commit ( &point_t );

	/* Custom operation for MPI_Reduce */
	MPI_Op op;
	MPI_Op_create( (MPI_User_function *)combineMeans, 1, &op );

	/* Print output to file */
	if ROOT
		fprintf ( fptr,"Number of processes: %d\n", nProc );

	/* For computing the elbow point using Sum of Squared Error */
	double oldDist, newDist, ratio, sse;
	int flag;

	/* Various timers */
	double totalTime[tStepTotal], proTime[tStepTotal], preProTime[tStepTotal];

	for ( i = 0; i < tStepTotal; i++ ) {
		totalTime[i] = 0;
		proTime[i] = 0;
		preProTime[i] = 0;
	}

	/* ========== Distribute tStepTotal time steps among nProc processes ========== */
	
	tpProc = tStepTotal / nProc;
	tOffset = myRank * tpProc;

	/* Ex: Distributes 11 time steps into 3, 4, 4 time steps among 3 processes */
	int val = nProc - (tStepTotal % nProc);
	if (tStepTotal % nProc != 0) {
		if ( myRank > val - 1 ) {
			tOffset = val * tpProc + ( myRank - val ) * (tpProc + 1 );
			tpProc++;
		}
	}

	#if DEBUG == 0
		printf ( "[%d] tpProc = %d, tOffset = %d\n", myRank, tpProc, tOffset );
	#endif

	/* For Collectives (output generation) */
	int kNosSendBuf[tpProc], kNosRecvBuf[tStepTotal], kNosRecvCounts[nProc], kNosDisPls[nProc], kValsDisPls[nProc], step = 0;
	index = KMAX * tpProc;
	point kValsSendBuf[index];

	/* Cheaper than having to allocate dynamic memory every single time in each of the multiple iterations */
	point means[KMAX], pMeans[KMAX], oldMeans[KMAX];
	point *kValsRecvBuf;
	meanInd = 0;

	if ROOT {
		for ( i = 0; i < val; i++ )
			kNosRecvCounts[i] = tpProc;
		
		if (tStepTotal % nProc != 0)
			for ( i = val; i < nProc; i++ )
				kNosRecvCounts[i] = tpProc + 1;
		else
			for ( i = val; i < nProc; i++ )
				kNosRecvCounts[i] = tpProc;

		kNosDisPls[0] = 0;
		for ( i = 0; i < nProc; i++ )
			if ( i < nProc - 1 )
				kNosDisPls[i + 1] = kNosDisPls[i] + kNosRecvCounts[i];
	}

	#if DEBUG == 0
		printf ( "[%d] tpProc = %d, tOffset = %d, tStep[0] = %d\n", myRank, tpProc, tOffset, tStep );
	#endif

	/* Common across time steps */
	_totalTime = (MPI_Wtime() - _totalTime) / tpProc;

	/* Common across time steps */
	_timePP = ( MPI_Wtime() - _timePP ) / tpProc;
	
	sse = 0;

	for ( tStep = tOffset; tStep < tOffset + tpProc; tStep++ ) {
		
		/* Start timers for this time step */
		totalTime[tStep] = MPI_Wtime();
		preProTime[tStep] = MPI_Wtime();

		/* Prepare actual data file path */
		sprintf ( path, "%s/file%d", dPath, tStep );

		#if DEBUG == 0
			if ROOT
				printf("[%d] %s\n", myRank, path);
		#endif

		/* Open file for read */
		error = MPI_File_open ( MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &dFile );
		if (error) {
			if ROOT
				printf ("Error reading from file\n");
			MPI_Finalize();
			exit(1);
		}

		#if DEBUG == 0
			if ROOT
				printf ( "[%d] Read file%s succesfully \n", myRank, dPath);
		#endif

		/* File (data) Size in bytes */
		MPI_File_get_size ( dFile, &dfSize );

		/* Every point is represented by 4 doubles */
		nPoints = dfSize / ( 4 * DOUBLE );

		pdBuf = malloc ( dfSize );
		pPoints = malloc ( POINT * nPoints );
		if ROOT {
			if ( pdBuf == NULL || pPoints == NULL)
				printf ("Error in memory allocation\n");
		}

		/* Read file data into buffer */
		MPI_File_read ( dFile, pdBuf, nPoints * 4, MPI_DOUBLE, MPI_STATUS_IGNORE );

		#if DEBUG == 0
			if ROOT
				printf ( "MPI_File_closed on %s\n", path);
		#endif

		MPI_File_close ( &dFile );

		/* Read all points */
		for ( i = 0; i < nPoints; i++ ) {
			pPoints[i].info = pdBuf[i * 4];
			pPoints[i].x = pdBuf[i * 4 + 1];
			pPoints[i].y = pdBuf[i * 4 + 2];
			pPoints[i].z = pdBuf[i * 4 + 3];
		}

		preProTime[tStep] = (MPI_Wtime() - preProTime[tStep]) + _timePP;

		oldDist = INT_MAX;
		flag = 0;

		proTime[tStep] = MPI_Wtime();

		/* K-Means */

		if ( tStepTotal == D1 )
			k = d1kSteps[tStep];
		else
			k = d2kSteps[tStep];

		#if DEBUG == 0
			if ROOT
				printf ( "[%d] Running for k = %d\n", myRank, k);
		#endif

		threshold = 100;

		/* Set seed for pseudo random no generation - input current time for making it random */
		srand(time(NULL));

		/* Initialize random cluster means *
		for( i = 0; i < k; i++) {
			index = rand() % nPoints;
			means[i].x = pPoints[index].x;
			means[i].y = pPoints[index].y;
			means[i].z = pPoints[index].z;
			means[i].info = 0;
		}*/

		/* Initialize Cluster Means */
		kMeansPP ( pPoints, nPoints, means, k);

		#if DEBUG == 0
			if ROOT
				print ( "Initial means:", k, means);
			//printf ( "Starting Parallel k-means %d\n", myRank );
		#endif

		iters = 0;

		/* Iterations of k-means */
		while ( threshold > 0.0001 && iters < 200 ) {
			iters++;
			
			#if DEBUG == 0
				printf ( "Done MPI_BCast %d\n", myRank );
				printf ("%d ", myRank);
				print ("Received means:", k, means );
			#endif

			zeroMeans ( k, pMeans );

			newDist = 0;
			/* Assign points to clusters */
			for ( i = 0; i < nPoints; i++ ) {
				newDist += clusterPoint ( &pPoints[i], k, means, pMeans );
			}

			#if DEBUG == 0
				printf("%d ", myRank);
				print ("Computed means:", k, pMeans );
			#endif

			if ROOT {
				for ( i = 0; i < k; i++ ) {
					oldMeans[i].x = means[i].x;
					oldMeans[i].y = means[i].y;
					oldMeans[i].z = means[i].z;
					oldMeans[i].info = means[i].info;
				}
				zeroMeans (k, means );
			}

			for ( i = 0; i < k; i++ ) {
				n = pMeans[i].info;
				if ( n != 0 ) {
					means[i].x = pMeans[i].x / n;
					means[i].y = pMeans[i].y / n;
					means[i].z = pMeans[i].z / n;
					means[i].info = n;
				}
			}
			
			#if DEBUG == 0
				if (myRank == 0)
					print ("Reduced means:", k, means );
			#endif

			if ROOT {
				threshold = 0;
				for ( i = 0; i < k; i++ )
					threshold += distSq ( oldMeans[i], means[i] );
			}
		}
				
		proTime[tStep] = (MPI_Wtime() - proTime[tStep]);
		
		#if DEBUG == 0
			if ROOT
				print ("New means:", k, means );
		#endif

		/* Disabled optimal k finding
		ratio = newDist / oldDist;
		if ( ratio > 0.95 ) {
			if ( flag == 0 )
				break;
			else
				flag--;
		}
		sse = 0;
		for ( i = 0; i < nPoints; i++ ) {
			index = (int)(pPoints[i].info);
			sse += distSq ( means[index], pPoints[i] );
		}
		*/
		
		oldDist = newDist;
		newDist = 0;

		if ( k > KMAX )
			k = KMAX;
		
		int kPrime = 0;

		for ( i = 0 ; i < k; i++ ) {
			/* Filter out bad means with zero membership */
			if ( means[i].info != 0 ) {
				kValsSendBuf[meanInd++] = means[i];
				kPrime++;
			}
		}

		k = kPrime;

		kNosSendBuf[step] = k;
		step++;

		totalTime[tStep] = (MPI_Wtime() - totalTime[tStep]) + _totalTime;
	}

	/* Gather # clusters for each of tSteps */
	MPI_Gatherv ( kNosSendBuf, tpProc, MPI_INT, kNosRecvBuf, kNosRecvCounts, kNosDisPls, MPI_INT, 0, MPI_COMM_WORLD );

	int sumkNosRecvCounts[nProc];
	int disp;
	int j;
	
	if ROOT {
		/* Calculate total # of cluster points (across time steps) per process */
		for ( i = 0; i < nProc; i++ ) {
			sumkNosRecvCounts[i] = 0;
			disp = kNosDisPls[i];
			for ( j = 0; j < kNosRecvCounts[i]; j++ )
				sumkNosRecvCounts[i] += kNosRecvBuf[disp + j];
		}

		/* Based on received values prepare buffers and displacements */
		kValsDisPls[0] = 0;
		for ( i = 0; i < nProc; i++ )
			if ( i < nProc - 1 )
				kValsDisPls[i + 1] = kValsDisPls[i] + sumkNosRecvCounts[i];

		#if DEBUG == 0
			printf ("\n");
			for ( i = 0; i < nProc; i++ ) {
				printf (" %d",sumkNosRecvCounts[i]);
			}
			printf ("\n");
			
			printf ("\n");
			for ( i = 0; i < nProc; i++ ) {
				printf (" %d",kValsDisPls[i]);
			}
			printf ("\n");	
		#endif

		int total = 0;
		for ( i = 0; i <  nProc; i++ ) {
			total += sumkNosRecvCounts[i];
		}

		kValsRecvBuf = malloc (POINT * total);	
	}

	/* kValsRecvBuf contains all of the means across times steps */
	MPI_Gatherv ( kValsSendBuf, meanInd, point_t, kValsRecvBuf, sumkNosRecvCounts, kValsDisPls, point_t, 0, MPI_COMM_WORLD );
	
	int tCount;
	meanInd = 0;
	if ROOT {
		/* Write results to file */	
		/* Scan kNosRecvBuf to get mean counts in each time step */
		for ( i = 0; i < tStepTotal; i++ ) {
			tCount = kNosRecvBuf[i];
			fprintf ( fptr, "T%d: %d", i + 1, tCount);
			//print means of step i to file
			for ( j = 0; j < tCount; j++ ) {
				fprintf ( fptr, ", <%d, %lf %lf %lf>", (int)kValsRecvBuf[meanInd].info, kValsRecvBuf[meanInd].x, kValsRecvBuf[meanInd].y, kValsRecvBuf[meanInd].z);
				meanInd++;
			}
			fprintf ( fptr, "\n" );
		}
	}

	#if DEBUG == 0
		if ROOT {
			for ( i = 0; i < tStepTotal; i++ )
				printf ("t%d: %d\n", i, kNosRecvBuf[i] );
		}
	#endif

	if ROOT {
		MPI_Reduce ( MPI_IN_PLACE, preProTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce ( MPI_IN_PLACE, proTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce ( MPI_IN_PLACE, totalTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	}
	else {
		MPI_Reduce ( preProTime, preProTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce ( proTime, proTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce ( totalTime, totalTime, tStepTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	}

	if ROOT {
		#if DEBUG == 0
			printf ("Time_Step\tPre_Proc\tProc\t\tTotal\n");
			for ( i = 0; i < tStepTotal; i++ ) {
				printf ("%d\t\t%lf\t%lf\t%lf\n", i, preProTime[i], proTime[i], totalTime[i] );
			}
		#endif

		avgPPTime = avgTime ( kNosRecvCounts, nProc, preProTime, tStepTotal);
		avgPTime = avgTime ( kNosRecvCounts, nProc, proTime, tStepTotal);
		pTotalTime = nTotalTime( kNosRecvCounts, nProc, totalTime, tStepTotal);
		
		fprintf ( pptFptr, "%d,%lf\n", nProc, avgPPTime );
		fprintf ( ptFptr, "%d,%lf\n", nProc, avgPTime );
		fprintf ( ttFptr, "%d,%lf\n", nProc, pTotalTime );
		
		fprintf ( fptr, "\nAverage time to pre-process: %lf\n", avgPPTime );
		fprintf ( fptr, "Average time to process: %lf\n", avgPTime );
		fprintf ( fptr, "Total time: %lf\n", pTotalTime );
	}

	if ROOT {
		fclose (fptr);
		fclose (pptFptr);
		fclose (ptFptr);
		fclose (ttFptr);
	}

	MPI_Type_free (&point_t);
	MPI_Op_free (&op);

	MPI_Finalize();
}