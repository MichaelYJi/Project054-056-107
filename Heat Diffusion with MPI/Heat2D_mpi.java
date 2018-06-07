//******************************************************************************
//**  File: Heat2D_mpi.java
//**  Modified by: Michael Ji
//******************************************************************************

import mpi.*;
import java.util.Date;
import java.util.ArrayList;

public class Heat2D_mpi
{
	private static double a = 1.0;  // heat speed
	private static double dt = 1.0; // time quantum
	private static double dd = 2.0; // change in

	public static void main(String[] args) throws MPIException
	{
    	MPI.Init(args);

    	int size = Integer.parseInt( args[1]);			// size * size
    	int max_time = Integer.parseInt( args[2]);		// # of cycles running
    	int heat_time = Integer.parseInt( args[3]);	// heated for this many cycles
    	int interval  = Integer.parseInt( args[4]);	// prints out intermediate status every # cycles
    	double r = a * dt / (dd * dd);

    	int mpiSize = MPI.COMM_WORLD.Size();

    	int tag = 0;

// provided code----------------------------------------------------------------
		// create a space
		double[][][] z = new double[2][size][size];
    	for ( int p = 0; p < 2; p++ )
    		for ( int x = 0; x < size; x++ )
    			for ( int y = 0; y < size; y++ )
    			{
					z[p][x][y] = 0.0; // no heat or cold
				}

		// start a timer
		Date startTime = new Date();

		// simulate heat diffusion
		for (int t = 0; t < max_time; t++)
		{
	    	int p = t % 2; // p = 0 or 1: indicates the phase

	    	// two left-most and two right-most columns are identical
	    	for (int y = 0; y < size; y++)
	    	{
	    		z[p][0][y] = z[p][1][y];
	    		z[p][size - 1][y] = z[p][size - 2][y];
	    	}

	    	// two upper and lower rows are identical
	    	for (int x = 0; x < size; x++)
	    	{
	    		z[p][x][0] = z[p][x][1];
	    		z[p][x][size - 1] = z[p][x][size - 2];
	    	}

	    	// keep heating the bottom until t < heat_time
	    	if (t < heat_time )
	   		{
	   			for ( int x = size /3; x < size / 3 * 2; x++ )
				z[p][x][0] = 19.0; // heat
			}
//end provided code-------------------------------------------------------------

			// boundary exchange
			// if more than 1 node, trade edge rows with stripe next door
			if (mpiSize > 1)
			{
				double arr[] = new double[size];
				double arr1[] = new double[size];

				// rank 0, left edge case
				if (MPI.COMM_WORLD.Rank() == 0)
				{
					// package up row to send
					int rightB = size / mpiSize;
					for (int y = 0; y < size; y++)
					{							
						arr[y] = z[p][rightB - 1][y];
					}

					// send right row
					MPI.COMM_WORLD.Send(arr, 0, size, MPI.DOUBLE, 1, tag);

					// receive right row
					MPI.COMM_WORLD.Recv(arr1, 0, size, MPI.DOUBLE, 1, tag);

					// implement into right row
					for (int y = 0; y < size; y++)
					{
						z[p][rightB][y] = arr1[y];						
					}
				}

				// last rank, right edge case
				else if (MPI.COMM_WORLD.Rank() == mpiSize - 1)
				{
					// boundary calculation
					double rawBoundary = (double) size / (double) mpiSize * (double) MPI.COMM_WORLD.Rank();
					int boundary = (int) rawBoundary;

					// receive left row
					MPI.COMM_WORLD.Recv(arr1, 0, size, MPI.DOUBLE, mpiSize - 2, tag);

					// implement into right row
					for (int y = 0; y < size; y++)
					{
						z[p][boundary - 1][y] = arr1[y];
					}

					// package up row
					for (int y = 0; y < size; y++)
					{
						arr[y] = z[p][boundary][y];
					}
					// send left row
					MPI.COMM_WORLD.Send(arr, 0, size, MPI.DOUBLE, mpiSize - 2, tag);					
				}
				else // inbetween cases
				{
					// left boundary calculation
					double rawLeftBoundary = (double) size / (double) mpiSize * (double) MPI.COMM_WORLD.Rank();
					int leftBoundary = (int) rawLeftBoundary;
					
					// package up left row
					for (int y = 0; y < size; y++)
					{
						arr[y] = z[p][leftBoundary][y];
					}

					// send left row
					MPI.COMM_WORLD.Send(arr, 0, size, MPI.DOUBLE, MPI.COMM_WORLD.Rank() - 1, tag);

					// receive left row
					MPI.COMM_WORLD.Recv(arr1, 0, size, MPI.DOUBLE, MPI.COMM_WORLD.Rank() - 1, tag);
					// implement into left row
					for (int y = 0; y < size; y++)
					{
						z[p][leftBoundary - 1][y] = arr1[y];
					}

					// right boundry calculation
					double rawRightBoundary = (double) size / (double) mpiSize * ((double) MPI.COMM_WORLD.Rank() + 1.0);
					int rightBoundary = (int) rawRightBoundary;

					// package up rght rowsend right row
					for (int y = 0; y < size; y++)
					{
						arr[y] = z[p][rightBoundary - 1][y];
					}

					// send right row
					MPI.COMM_WORLD.Send(arr, 0, size, MPI.DOUBLE, MPI.COMM_WORLD.Rank() + 1, tag);

					// receive right row
					MPI.COMM_WORLD.Recv(arr1, 0, size, MPI.DOUBLE, MPI.COMM_WORLD.Rank() + 1, tag);

					// implement into right row
					for (int y = 0; y < size; y++)
					{
						z[p][rightBoundary][y] = arr1[y];
					}
				}


			}
	
			// if time matches the interval, send stripes to rank 0
			if (interval != 0 && (t % interval == 0 || t == max_time - 1))
			{
				// if there is more than 1 node
				if (mpiSize > 1)
				{
					double rawAvg = size / mpiSize;

					// based of specs
					int stripeLength = 2 * size * size;
					double[] stripe = new double[stripeLength];

					// not rank 0, sends current stripe to rank 0
					if (MPI.COMM_WORLD.Rank() != 0)
					{
						// partition start & end
						int start = (int)(rawAvg * (double)MPI.COMM_WORLD.Rank());
						int end = (int)(rawAvg * ((double)MPI.COMM_WORLD.Rank() + 1.0));

						// converting the stripe, 3d array into 1d
						for (int x = start; x < end; x++)
						{
							for (int y = 0; y < size; y++)
							{
								stripe[p * size * size + x * size + y] = z[p][x][y];
							}
						}
						// calculating the length of the array
						int rangeStart = p * size * size + start * size + 0;
						int rangeEnd = (int)rawAvg * size;

						// sending sripe converted to 1d array to rank 0
						MPI.COMM_WORLD.Send(stripe, rangeStart, rangeEnd, MPI.DOUBLE, 0, tag);
					}

					//rank 0 will accept array data and will implement it
					if (MPI.COMM_WORLD.Rank() == 0)
					{
						// traverse through each rank
						for (int i = 1; i < mpiSize; i++)
						{
							// calculate the start and finish offset
							int start = (int)(rawAvg * (double)i);
							int end = (int)(rawAvg * ((double)i + 1.0));

							// calculating the length of the array
							int rangeStart = p * size * size + start * size + 0;
							int rangeEnd = (int)rawAvg * size;

							// receive array from each rank
							MPI.COMM_WORLD.Recv(stripe, rangeStart, rangeEnd, MPI.DOUBLE, i, tag);
							// convert each array back into 3d with equation from spec
							for (int x = start; x < end; x++)
							{
								for (int y = 0; y < size; y++)
								{
									z[p][x][y] = stripe[p * size * size + size * x + y];
								}
							}
						}
					}
				}
			}

// provided code----------------------------------------------------------------
			// only rank 0 needs to print this out
			// rank 0 must receive all stripes before reaching this
			// display intermediate results
			if (MPI.COMM_WORLD.Rank() == 0)
			{
				if (interval != 0 && (t % interval == 0 || t == max_time - 1))
				{
					System.out.println( "time = " + t );
					for (int y = 0; y < size; y++ )
					{
						for (int x = 0; x < size; x++ )
						{
							System.out.print((int)(Math.floor(z[p][x][y] / 2)) + " ");
						}
						System.out.println();
					}
					System.out.println();
				}
			}

			// there is only 1 node, heats
			if (mpiSize == 1)
			{
				int p2 = (p + 1) % 2;
				for (int x = 1; x < size - 1; x++)
				{
					for (int y = 1; y < size - 1; y++)
					{
						z[p2][x][y] = z[p][x][y] + r * (z[p][x + 1][y] - 2 * z[p][x][y] + z[p][x - 1][y]) + r * (z[p][x][y + 1] - 2 * z[p][x][y] + z[p][x][y - 1]);
					}
				}
			}
// end provided code------------------------------------------------------------

			// heating, perform forward Euler method if more than 1 node
			if (mpiSize > 1)
			{
				double avg = (double)size / (double)mpiSize;
				
				// partition start
				double rawLower = avg * (double)MPI.COMM_WORLD.Rank();
				int start = (int)rawLower;

				// partition end
				double rawUpper = avg * ((double)(MPI.COMM_WORLD.Rank() + 1.0));
				int end = (int)rawUpper;

				// heating for rank 0
				if (MPI.COMM_WORLD.Rank() == 0)
				{
					int p2 = (p + 1) % 2;
					for (int x = 1; x < end; x++) // heats from 1 to end of its partition
					{
						for (int y = 1; y < size - 1; y++)
						{
							z[p2][x][y] = z[p][x][y] + r * (z[p][x + 1][y] - 2 * z[p][x][y] + z[p][x - 1][y]) + r * (z[p][x][y + 1] - 2 * z[p][x][y] + z[p][x][y - 1]);
						}
					}
				}

				// heating for last rank
				else if (MPI.COMM_WORLD.Rank() == mpiSize - 1)
				{
					int p2 = (p + 1) % 2;
					for (int x = start; x < end - 1; x++) // heats from stripe start until end - 1
					{
						for (int y = 1; y < size - 1; y++)
						{
							z[p2][x][y] = z[p][x][y] + r * (z[p][x + 1][y] - 2 * z[p][x][y] + z[p][x - 1][y]) + r * (z[p][x][y + 1] - 2 * z[p][x][y] + z[p][x][y - 1]);
						}
					}
				}

				// heating for in between cases
				else
				{
					int p2 = (p + 1) % 2;
					for (int x = start; x < end; x++) // heats stripe start to end
					{
						for (int y = 1; y < size - 1; y++)
						{
							z[p2][x][y] = z[p][x][y] + r * (z[p][x + 1][y] - 2 * z[p][x][y] + z[p][x - 1][y]) + r * (z[p][x][y + 1] - 2 * z[p][x][y] + z[p][x][y - 1]);
						}
					}
				}
			}
		} // end of simulation

		// rank 0 prints out results
		if (MPI.COMM_WORLD.Rank() == 0)
		{
			// finish the timer
			Date endTime = new Date();
			System.out.println("Elapsed time = " + (endTime.getTime() - startTime.getTime()));
		}

		// Terminate the MPI library.
		MPI.Finalize();
	}
}