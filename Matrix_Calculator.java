import java.util.Calendar;

/**
 * @Name: Pareekshit Ravi
 * @MacID: ravip2
 * @Student_Number: 1407109
 * @Description: Multiply N number of matrices and take the inverse of the final product matrix. The determinant of the 
 * 				 product matrix has been calculated iteratively for extra credit
 */
public class Matrix_Calculator {

	public static void main(String[] args) {
		//number of matrices
		int N = Integer.parseInt(args[0]);
		//count for the number of matrices
		int count = 0;
		//count for which argument to read from
		int argsCount = 1;
		//3d array to hold all the information
		double[][][] allMatrix = new double [N][][];
		//getting the dimensions of each individual matrix
		for (int i=1;i<N*2+1;i=i+2){
			//extracting the dimensions from the arguments
			int m = Integer.parseInt(args[i]);
			int n = Integer.parseInt(args[i+1]);
			//initializing the size of each individual matrix
			allMatrix[count] = new double[m][n];
			//next matrix
			count++;
			//increment 2 for argsCount because read 2 arguments
			argsCount+=2;
		}
		//setting all the matrices' values
		//loop for number of matrices
		for(int i=0;i<N;i++){
			//getting the number of rows and columns based on the size of the array made earlier
			int col = allMatrix[i].length;
			int row = allMatrix[i][0].length;
			//looping through the i'th matrix for the rows and columns
			for(int j=0;j<col;j++){
				for(int k=0;k<row;k++){
					//setting the value of the i'th matrix at column m and row n
					allMatrix[i][j][k]=Integer.parseInt(args[argsCount]);
					//incrementing the counter because read 1 argument was read
					argsCount++;
				}
			}
		}
		//			printMatrix(allMatrix);
		//multiplying all the matrices with each other
		//the final matrix at index N-1, will have a matrix that is the product of all the previous matrices
		//doing matrix multiplication from the right
		for (int i=1;i<N;i++){
			//multiplying 2 matrices together, at index i and the one before it
			//also ensures that the following matrices are multiplied to the right
			double[][] M = multiplyMatrix(allMatrix,i);
			if (M==null){
				System.exit(0);
			}
			//redefining the size of the matrix after the multiplication is complete
			allMatrix[i] = new double[M.length][M[0].length];
			//saving the product matrix to the allMatrix variable so that is used for the following multiplications
			allMatrix[i] = M;
		}
		//identifying the final product matrix as the N-1 index
		double[][] productMatrix = allMatrix[N-1];
		//matrix made for debugging purposes
		//		double[][] productMatrix;
		//		//setting matrix for me to debug
		//		double[][] test = {
		//				{5,-4,5,0},
		//				{1,0,8,2},
		//				{0,7,-1,-4},
		//				{2,7,4,2},
		//		};
		//		productMatrix = test;
		//		System.out.println("-----------------------------------Product Matrix-----------------------------");
		//		printMatrix(productMatrix);
		//		System.out.println();
		//checking that the matrix is a square
		if (productMatrix.length!=productMatrix[0].length){
			//matrix cannot be inverted as it is not square
			System.out.println("Matrix not invertible");
		}
		//matrix is a square matrix
		else{
			//initializing the variable to hold the cofactor matrix
			double [][] cofactorMatrix = new double[productMatrix.length][productMatrix.length];
			//iterating through every term in the product matrix through the columns and rows
			for(int i=0;i<productMatrix.length;i++){
				for(int j=0;j<productMatrix.length;j++){
					//calculating the minor of each element in the product matrix
					double[][] minor = getMinorMatrix(productMatrix,i,j);
					//calculating the determinant of the minor matrix
					//the 1 that is passed is the value that the minor needs to be multiplied by, it does not affect the determinant here
					//the i+j+2 is the formula used to calculate the sign of the term, + or -
					double determinant = getDeterminant(minor,1,i+j+2);
					//setting the determinant in the cofactor matrix
					cofactorMatrix[i][j] = determinant;
				}
			}
			//			System.out.println("-----------------------------------Cofactor Matrix-----------------------------");
			//						printMatrix(cofactorMatrix);
			//						System.out.println();
			//transposing the cofactor matrix
			double [][] transposedMatrix = transposeMatrix(cofactorMatrix);
			//			System.out.println("-----------------------------------Transpose Matrix-----------------------------");
			//			printMatrix(transposedMatrix);
			//			System.out.println();
			//timing the difference between finding the determinant iteratively and recursively
			//the start time of finding iterative determinant
			long startTime = System.nanoTime();
			//method to find the determinant iteratively of the product matrix
			double matrixDet = getDeterminantIterative(productMatrix);
			//the matrix has been calculated and the time is recorded
			long endTime = System.nanoTime();
			//printing the time it took to calculate the determinant and printing it
			//			System.out.println("Iterative determinant: " + (long)((endTime-startTime)) + " nanoseconds");
			startTime = System.nanoTime();
			//calculating the determinant of the product matrix
			//the 1 is the value that the overall matrix needs to be multiplied by
			//the 0 because the sign is irrelevant for this calculation
			//			matrixDet = getDeterminant(productMatrix,1,0);
			endTime = System.nanoTime();
			//			System.out.println("Recursive determinant: " + (long)((endTime-startTime)) + " nanoseconds");
			//System.out.println("Determinant of product matrix is: " + matrixDet);
			//if the determinant is not 0, it is invertible
			if (matrixDet!=0 && !Double.isNaN(matrixDet)){
				//initializing the inverse matrix variable and the size to be the same as the transposed matrix
				double[][] inverseMatrix = new double [transposedMatrix.length][transposedMatrix.length];
				//iterating through the full transposed matrix to every element
				for(int i=0;i<transposedMatrix.length;i++){ //columns
					for(int j=0;j<transposedMatrix.length;j++){ //rows
						//dividing each element of the transposed matrix by the product matrix's determinant
						inverseMatrix[i][j] = transposedMatrix[i][j]/matrixDet;
					}
				}
				//				printMatrix(inverseMatrix);

				//printing the answer in a straight line
				printAnswer(inverseMatrix);
			}
			//det is 0 or NaN
			else{
				System.out.println("Matrix not invertible");
			}
		}
	}
	private static double[][] transposeMatrix(double[][] cofactorMatrix) {
		//initializing a matrix that will hold the transposed matrix that is the same size as the cofactor matrix
		double[][] temp = new double[cofactorMatrix.length][cofactorMatrix.length];
		//iterating through every term in the cofactor matrix
		for (int i = 0; i < cofactorMatrix.length; i++){ //columns
			for (int j = 0; j < cofactorMatrix.length; j++){ //rows
				//placing the term of the column in the row of the temp matrix to transpose it
				temp[j][i] = cofactorMatrix[i][j];
			}
		}
		return temp;
	}
	private static double[][] getMinorMatrix(double[][] productMatrix,int a,int b) {
		//initializing 2d array to hold the minor
		//the size is one less in the column and rows
		double[][] c = new double[productMatrix.length-1][productMatrix.length-1];
		//counters to hold the indexes of array c
		int x = 0; //rows
		int y = 0; //columns
		//if only a 1x1 matrix, return a 1d array with value 1, so when dividing by determinant later, the inverse is 1/det
		if (productMatrix.length==1){
			double[][] ans = {{1}};
			c = ans;
			return c;
		}
		//iterating through all the terms in the product matrix columns and rows
		for(int i=0;i<productMatrix.length;i++){
			for(int j=0;j<productMatrix.length;j++){
				//if the row and column indices are not the same as the one passed
				//both must be different in order for it to be chosen as part of the minor
				if (i!=a && j!=b){
					//populating the minor
					c [x][y] = productMatrix[i][j];
					//next element in the minor
					y++;
					//the row is done, next row
					if (y==(productMatrix.length-1)){
						//columns set back to first
						y=0;
						//next row
						x++;
					}
				}
			}
		}
		return c;
	}
	//a recursive function that constantly gets the minor of a matrix till it reaches a 2x2 matrix
	private static double getDeterminant(double[][] square, double x, int exponent){
		//initializing the determinant
		double det = 0;
		//case if there is a 1x1 matrix
		if(square.length==1){
			//the determinant is the value in the matrix itself
			det = square[0][0];
		}
		//main base case for recursion
		else if(square.length==2){
			//calculating the determinant of a 2x2 matrix
			det = (square[0][0]*square[1][1])-(square[0][1]*square[1][0]);
		}
		else{
			//iterating through the first row of the matrix given
			for(int i =0;i<square.length;i++){
				//calculating the minor of that element
				//the minor will be 1 less size in length and width
				double[][] temp = getMinorMatrix(square,i,0);
				//getting the determinant of the smaller matrix and taking the sum of all the determinants of the minor matrices from that row
				det+=getDeterminant(temp,square[i][0], i+2);
			}
		}
		//when it has finished calculating the determinant of any sized matrix, it multiplies that by the coefficient that the matrix was a minor of
		det*=x;
		//calculating if the sign for this term is positive or negative with this formula
		int sign = (int) Math.pow(-1,exponent);
		//setting the sign
		det *=sign;
		return det;
	}
	private static double[][] multiplyMatrix(double[][][] allMatrix, int i) {
		//Initializing a 2d array for the product
		double[][] temp = new double[allMatrix[i-1].length][allMatrix[i][0].length];
		//		System.out.println(allMatrix[i-1][0].length + " " + allMatrix[i].length);
		//making sure it can be multiplied
		if (allMatrix[i-1][0].length==allMatrix[i].length){
			//columns
			for(int j=0;j<temp.length;j++){
				//row
				for(int k=0;k<temp[0].length;k++){
					//executing the matrix multiplication by multiplying every element in a row of matrix b by a column in matrix a
					for(int l=0;l<allMatrix[i].length;l++){
						//iterating through the two arrays and multiplying the terms
						//taking the sum of the product of the terms from the row * column of the 2 matrices
						temp[j][k]+= allMatrix[i-1][j][l] * allMatrix[i][l][k];
					}
				}
			}
			//			printMatrix(temp);
			//product is calculated and returned
			return temp;
		}
		//the columns of the left matrix is not the same as the rows as the right matrix thus it cannot be multiplied
		else{
			System.out.println("Multiplication Error");
			return null;
		}
	}
	private static double getDeterminantIterative(double[][] A) {
		//number of columns
		int N  = A.length;
		//tracking the row swaps
		int globalFactor = 1;
		//iterating through the matrix from the rows
		for (int a = 0; a < N; a++) {
			// find pivot column and swap
			//the max value in the row
			int max = a;
			//checking the other values in the same row
			for (int i = a+1; i < N; i++) {
				//if the value in another row is larger
				if (Math.abs(A[i][a]) > Math.abs(A[max][a])) {
					//setting the column with the max value in a row
					max = i;
				}
			}
			//swapping the columns
			double[] temp = A[a]; 
			A[a] = A[max]; 
			A[max] = temp;
			//if the column does swap. Max was initially set to a, but if after the for loop above, it never changed, it still swaps above, thus
			//if the values are different, it means that they really did swap and not that it was the same column replaced with the same column
			if (a!=max){
				//notifying that it has swapped columns, thus the determinant needs to by changed sign
				globalFactor*=-1;
			}
			//iterating through the columns
			for (int i = a+1; i < N; i++) {
				//the factor to multiply by to zero other values of column
				double factor = A[i][a] / A[a][a];
				// zeroing the other elements of the column below making a trianglar matrix
				for (int j = a; j < N; j++) {
					A[i][j] -= factor * A[a][j];
				}
			}
		}
		//		printMatrix(A);
		//initializing the determinant value
		double det = 1;
		//iterating through the triangular matrix for the columns and the rows
		for (int i=0;i<A.length;i++){
			for (int j=0;j<A[i].length;j++){
				//at the diagonal
				if (i==j){
					//an error for the value, if it is below 0.00001, then just treat it as a 0
					if ((Math.abs(A[i][j])<1e-5)){
						det*=0;
					}
					//at the diagonal of the matrix, product of the diagonals is the determinant
					det*=A[i][j];
				}
			}
		}
		//divide by the global factor because of the number of column swaps
		det/=globalFactor;
		//if the det is small enough, make it 0 similar to the error check above
		if ((Math.abs(det)<1e-5 && Math.abs(det)>0)){
			return 0;
		}
		else{
			return det;
		}
	}

	@SuppressWarnings("unused")
	private static void printMatrix(double[][][] allMatrix) {
		//iterating through all the elements of the array
		//going through all the matrices
		for(int i=0;i<allMatrix.length;i++){
			//the columns fo the matrix
			for(int j=0;j<allMatrix[i].length;j++){
				//rows of the matrix
				for(int k=0;k<allMatrix[i][j].length;k++){
					//printing each element of the matrix with a space
					System.out.print(allMatrix[i][j][k] + " ");
				}
				//going to the next line because the row is done
				System.out.println();
			}
			//going 2 lines down to represent the completion of that matrix
			System.out.println("/n");
		}
	}
	@SuppressWarnings("unused")
	private static void printMatrix(double[][] productMatrix) {
		//iterating through the full matrix
		for(int j=0;j<productMatrix.length;j++){
			for(int k=0;k<productMatrix[j].length;k++){
				//printing the element and a space
				System.out.print(productMatrix[j][k] + " ");
			}
			//printing a line to go to the next row
			System.out.println();
		}
	} 
	private static void printAnswer(double[][] inverseMatrix) {
		//iterating through the full inverse matrix columns and rows
		for(int j=0;j<inverseMatrix.length;j++){
			for(int k=0;k<inverseMatrix[j].length;k++){
				//printing only 2 decimal places of the matrix
				System.out.printf("%.2f",inverseMatrix[j][k]);
				//adding a space between each value
				System.out.print(" ");
			}
		}
	}
}