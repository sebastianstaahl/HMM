import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.lang.Math;
import java.util.Collections;


public class HMM2{

public static double [][] matrixMultiplication(double [][] Matrix1,double [][] Matrix2 ){
	int r1 = Matrix1.length;
	int c1 = Matrix1[0].length;
	int r2 = Matrix2.length;
	int c2 = Matrix2[0].length;

	double[][] Matrixres = new double [r1][c2];

	for (int i = 0; i < r1; i++) {
		for(int j = 0; j < c2; j++) {
		 for(int k = 0; k < c1; k++) {
			Matrixres[i][j] += Matrix1[i][k] * Matrix2[k][j];
		 }
		}
	}return Matrixres;

}



public static int[]  observationFiller(String inputLine){
	String [] lineArray = inputLine.split("\\s+");  //text.split("\\s+");

	//String[] cloneOfArray = Arrays.copyOf(lineArray, lineArray.length);
	//System.out.println(Arrays.toString(cloneOfArray));

	int noelements,itterator;
	noelements=Integer.parseInt(lineArray[0]);

	itterator=1;

	int[] vectorF=new int [noelements];

	//int i = startingNumber; i <= 100; i++
	for(int i = 0; i < vectorF.length; i++) //Denna gar igenom varje rad och for varje rad vill vi ga igen varje kolumn
	{
				vectorF[i]=Integer.parseInt(lineArray[itterator]);
				itterator++;

	}
	//String[] newArray = Arrays.copyOfRange(line, startIndex, endIndex);
	return vectorF;

}




public static double[][]  matrixFiller(String inputLine){
	String [] lineArray = inputLine.split("\\s+");  //text.split("\\s+");

	String[] cloneOfArray = Arrays.copyOf(lineArray, lineArray.length);
	//System.out.println(Arrays.toString(cloneOfArray));

	int rows,cols,itterator;
	rows=Integer.parseInt(lineArray[0]);
	cols=Integer.parseInt(lineArray[1]);
	//System.out.println(rows);
	//System.out.println(cols);

	itterator=2;

	double[][] matrixF=new double [rows][cols];

	//int i = startingNumber; i <= 100; i++
	for(int r = 0; r < matrixF.length; r++) //Denna gar igenom varje rad och for varje rad vill vi ga igen varje kolumn
	{
		//System.out.println("current row:"+r);
		for(int c = 0; c < matrixF[r].length; c++) // har gar vi igenom varje element i en rad.
		{
				//System.out.println("current col:"+c);
				//System.out.println(lineArray[itterator]);
				matrixF[r][c]=Double.parseDouble(lineArray[itterator]);
				itterator++;

		}
	}
	//String[] newArray = Arrays.copyOfRange(line, startIndex, endIndex);
	return matrixF;

}


public static double maxFinder(double [] vector1){
double biggest=vector1[0];
for(int i=1; i<vector1.length;i++){
	if(vector1[i]>biggest){
		biggest=vector1[i];
	}
} return biggest;
}

public static double [] getColumn(double [][] matrix1,int c){
	double [] extractVector =new double[matrix1.length];
	for (int i=0; i<matrix1.length;i++){
		extractVector[i]=matrix1[i][c];

	}return extractVector;
}


public static int indexOfMaxFinder(double [] vector1){
double biggest=vector1[0];
int indexOfBiggest=0;
for(int i=1; i<vector1.length;i++){
	if(vector1[i]>biggest){
		biggest=vector1[i];
		indexOfBiggest=i;
	}
} return indexOfBiggest;
}


public static void main (String[] args){

				String fileName = args[0];
				File file = new File(fileName);
				Scanner input;

				try{
					input = new Scanner(file);
					//while(input.hasNextLine()){
					//int row, col;
					String line=input.nextLine();
					double [][] transitionMatrix =matrixFiller(line);

					//int[] cloneOfArray = Arrays.copyOf(transitionMatrix[2], transitionMatrix[2].length);
					//System.out.println(Arrays.toString(cloneOfArray));
					//System.out.println("Transistionmatrix");


					String line2=input.nextLine();
					double [][] emissionMatrix = matrixFiller(line2);


					String line3=input.nextLine();
					double [] initialStateProb = matrixFiller(line3)[0];


					String line4=input.nextLine();
					int [] observationSequence = observationFiller(line4);

					//System.out.println(Arrays.toString(observationSequence));
					int firstobsvalandindex= observationSequence[0];

					int noStates = transitionMatrix[0].length;
					//System.out.println("nostates: "+noStates);
					int noTimeSteps =observationSequence.length;

					double []delta0=new double [noStates];
					int [][]indiciesMaxProb=new int [noStates][noTimeSteps];

					for(int i=0;i<noStates;i++){
						delta0[i]=Math.log(initialStateProb[i]*emissionMatrix[i][firstobsvalandindex]);
						indiciesMaxProb[i][0]=0;
					}


					//int indexMaxTimeStepProb=indexOfMaxFinder(delta0);
					//indiciesMaxProb[0]=indexMaxTimeStepProb;

					double [][]deltat=new double [noStates][noTimeSteps];

					for(int l=0; l<delta0.length;l++){
						deltat[l][0]=delta0[l];

					}

					for (int t =1; t <noTimeSteps;t++){
						//double [] tmpTimeStepProb = new double [noStates];
						//System.out.println("Timestep : "+t);
						for (int i=0;i<noStates;i++){
							double []tmpStateProb =new double[noStates];
							//System.out.println("i: "+i);
							for (int j=0;j<noStates;j++){
							//double maxDeltaPrev = Collections.max(deltat[t-1])
							//observationSequence[t] the value of this is the same as the index
							//System.out.println("i:"+i);
							//System.out.println("j:"+j);
							//System.out.println("deltat[j][t-1]"+deltat[j][t-1]);
							//System.out.println("Math.log(transitionMatrix[j][i])"+Math.log(transitionMatrix[j][i]));
							//System.out.println("Math.log(emissionMatrix[i][observationSequence[t]])"+Math.log(emissionMatrix[i][observationSequence[t]]));
							//System.out.println();
							tmpStateProb[j]=deltat[j][t-1]+Math.log(transitionMatrix[j][i])+Math.log(emissionMatrix[i][observationSequence[t]]);

							}
							//if(t==3){
							//System.out.println("State: "+i);
							//System.out.println(Arrays.toString(tmpStateProb));
							//System.out.println(maxFinder(tmpStateProb));
							deltat[i][t]=maxFinder(tmpStateProb);
							indiciesMaxProb[i][t]=indexOfMaxFinder(tmpStateProb);
							//System.out.println();
						}


					}
					double maxProbLastTimeStep = maxFinder(getColumn(deltat,noTimeSteps-1));
					int pointerLastTimeStep = indexOfMaxFinder(getColumn(deltat,noTimeSteps-1));

					//System.out.println(Arrays.deepToString(indiciesMaxProb));
					String res=Integer.toString(pointerLastTimeStep);
					//System.out.println(res);
					int [] resArray = new int [noTimeSteps];
					resArray[noTimeSteps-1]=pointerLastTimeStep;
					//System.out.println(Arrays.toString(resArray));

					for (int t =noTimeSteps-1; t> 0;t--){
						//System.out.println("t :"+t);
						resArray[t-1]=indiciesMaxProb[pointerLastTimeStep][t];
						pointerLastTimeStep=indiciesMaxProb[pointerLastTimeStep][t];
						//System.out.println(Arrays.toString(resArray));
						//System.out.println(res);
					}
					String val= Integer.toString(resArray[0]);

					for (int i=1; i<resArray.length;i++){
						val+=" "+Integer.toString(resArray[i]);
						//System.print(val)
					}
					System.out.println(val);
					//String resFinal = Arrays.toString(resArray);
					//System.out.println(String.join(resFinal));


				} catch (FileNotFoundException exception) {
					System.out.println("The file does not exist");}



		}

}
