import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;
import java.util.*;
import java.util.Arrays;
import java.lang.Math;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.DoubleStream;
import static java.lang.System.*;
import java.util.Random;


public class HMM{
  static final Random noint = new Random();
  // Good seeds: 397169450 816905934 1033159171 1116287974 32

  //1483962880 1448742097  789487785 1019674898

//noint.nextInt()

  static final int rint=Math.abs(noint.nextInt());//2069414876);//1657773347);//1019674898      100p: 1657773347


  static double epsilon =(Double) Math.pow(10, -50);
  static double []C;
  static double [][] transitionMatrix;
	static double [][] transitionMatrixCorrect;
  static double [][] emissionMatrix;
	static double [][] emissionMatrixCorrect;
  static double [] initialStateProb;
	static double [] initialStateProbCorrect;
  static double [][][] diGammat;
  static double [][] gammat;

  static int noStates = 3;
  static int noEmissions = 4;
  static double denom;
  static double numer;

  public HMM() {

		//String line=input.nextLine();
		transitionMatrix = new double[][] {{0.54,0.26,0.20}, {0.19,0.53,0.28}, {0.22,0.18,0.6}};
		transitionMatrixCorrect = new double [][] {{0.7, 0.05, 0.25}, {0.1, 0.8, 0.1}, {0.2, 0.3, 0.5}};


		//String line2=input.nextLine();
		 emissionMatrix = new double [][] {{0.5,0.2,0.11,0.19}, {0.22,0.28,0.23,0.27}, {0.19,0.21,0.15,0.45}};
		 emissionMatrixCorrect = new double[][]{{0.7, 0.2 , 0.1, 0}, {0.1, 0.4, 0.3, 0.2}, {0, 0.1, 0.2, 0.7}};


		//String line3=input.nextLine();
		 initialStateProb = new double [] {0.3,0.2,0.5};
		 initialStateProbCorrect = new double[]{1, 0, 0};

  }

	public static int[] observationFiller(String inputLine){
		String [] lineArray = inputLine.split("\\s+");  //text.split("\\s+");

		//String[] cloneOfArray = Arrays.copyOf(lineArray, lineArray.length);
		//System.out.println(Arrays.toString(cloneOfArray));

		int noelements,itterator;
		noelements=Integer.parseInt(lineArray[0]);

		itterator=1;

		int[] vectorF=new int [noelements];

		//int i = startingNumber; i <= 100; i++
		for(int i = 0; i < vectorF.length; i++) //Denna går igenom varje rad och för varje rad vill vi gå igen varje kolumn
		{
					vectorF[i]=Integer.parseInt(lineArray[itterator]);
					itterator++;

		}
		//String[] newArray = Arrays.copyOfRange(line, startIndex, endIndex);
		return vectorF;

	}


public static  double [][] alfaPassNO(int [] obs){
  double[] cloc = new double[obs.length];
  double [][] alfaloc = new double [noStates][obs.length];
  cloc[0]=0; //epsilon

  for (int i = 0; i < noStates; i++) {
    alfaloc[i][0]=initialStateProb[i]*emissionMatrix[i][obs[0]];
    cloc[0] += alfaloc[i][0];
  }
  cloc[0] = 1/(cloc[0]+epsilon);

  for (int i = 0; i < noStates; i++) {
    alfaloc[i][0]=cloc[0]*alfaloc[i][0];
  }

  for (int t = 1; t < obs.length; t++) {
    cloc[t]=0;
    for (int i = 0; i < noStates; i++) {
      cloc[i] = 0;
      for (int j = 0; j < noStates; j++) {
        alfaloc[i][t] =alfaloc[i][t]+ alfaloc[j][t-1]*transitionMatrix[j][i];
      }
      alfaloc[i][t]=alfaloc[i][t]*emissionMatrix[i][obs[t]];
      cloc[t] =cloc[t]+ alfaloc[i][t];
    }
    cloc[t] = 1/(cloc[t]+epsilon);

    for (int i = 0; i < noStates; i++) {
      alfaloc[i][t]=cloc[t]*alfaloc[i][t];
    }
  }
  return alfaloc;
}

        public static  double [][] matrixMultiplication(double [][] Matrix1,double [][] Matrix2 ){
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

public static  double [] getColumn(double [][] matrix1,int c){
    double [] extractVector =new double[matrix1.length];
    for (int i=0; i<matrix1.length;i++){
        extractVector[i]=matrix1[i][c];

    }return extractVector;
}

public static double alfaPassNOprob(int [] obs){
  double[] cloc = new double[obs.length];
  double [][] alfaloc = new double [noStates][obs.length];
  cloc[0]=0; //epsilon

  for (int i = 0; i < noStates; i++) {
    alfaloc[i][0]=initialStateProb[i]*emissionMatrix[i][obs[0]];
    cloc[0] += alfaloc[i][0];
  }
  cloc[0] = 1/(cloc[0]+epsilon);

  for (int i = 0; i < noStates; i++) {
    alfaloc[i][0]=cloc[0]*alfaloc[i][0];
  }

  for (int t = 1; t < obs.length; t++) {
    cloc[t]=0;
    for (int i = 0; i < noStates; i++) {
      cloc[i] = 0;
      for (int j = 0; j < noStates; j++) {
        alfaloc[i][t] =alfaloc[i][t]+ alfaloc[j][t-1]*transitionMatrix[j][i];
      }
      alfaloc[i][t]=alfaloc[i][t]*emissionMatrix[i][obs[t]];
      cloc[t] =cloc[t]+ alfaloc[i][t];
    }
    cloc[t] = 1/(cloc[t]+epsilon);

    for (int i = 0; i < noStates; i++) {
      alfaloc[i][t]=cloc[t]*alfaloc[i][t];
    }
  }
  double logProb = 0;

  for(int t=0; t<obs.length;t++){
    //if(cloc[t] > 0){
      logProb += Math.log(cloc[t]+epsilon);
  //}
}
logProb = -logProb;
  return logProb;
}

public static double[][] transposeMatrix(double [][] matrixToTrans){
    double[][] transposedMatrix = new double[matrixToTrans[0].length][matrixToTrans.length];
    for (int i = 0; i < matrixToTrans.length; i++)
        for (int j = 0; j < matrixToTrans[0].length; j++)
            transposedMatrix [j][i] = matrixToTrans[i][j];
    return transposedMatrix;
}



public static  double [][] alfapassScaling(int [] observationSequence,double [][] alfat){
            C[0]=0; //epsilon
                    for (int i = 0; i < noStates; i++) {
                        alfat[i][0]=initialStateProb[i]*emissionMatrix[i][observationSequence[0]];
              C[0] += alfat[i][0];
                    }
                    C[0] = 1/(C[0]+epsilon);

                    for (int i = 0; i < noStates; i++) {
                        alfat[i][0]=C[0]*alfat[i][0];
                    }

                    for (int t = 1; t < observationSequence.length; t++) {
                        C[t]=0;
                        for (int i = 0; i < noStates; i++) {
                            alfat[i][t] = 0;
                            for (int j = 0; j < noStates; j++) {
                                alfat[i][t] =alfat[i][t]+ alfat[j][t-1]*transitionMatrix[j][i];
                            }
                            alfat[i][t]=alfat[i][t]*emissionMatrix[i][observationSequence[t]];
                            C[t] =C[t]+ alfat[i][t];
                        }
                        C[t] = 1/(C[t]+epsilon);

                        for (int i = 0; i < noStates; i++) {
                            alfat[i][t]=C[t]*alfat[i][t];
                        }
                    }
                    return alfat;
                }


          public static  double [][] betapass( int [] observationSequence ,double [][]betat){

             for (int i = 0; i < noStates; i++) {
               betat[i][observationSequence.length-1] = C[observationSequence.length-1];
             }

             for (int t = observationSequence.length-2; t > -1 ; t--){
               for (int i = 0; i < noStates; i++) {
                 betat[i][t] = 0;
                 for (int j = 0; j < noStates; j++) {
                   betat[i][t] =betat[i][t]+ transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1];
             }
                 betat[i][t] = C[t]*betat[i][t];
               }
             }
             return betat;
           }

           public static void diGammaGamma (double [][] betat,double [][] alfat, int [] observationSequence){
             for (int t = 0; t < observationSequence.length-1; t++) {
               denom = 0;
               for (int i = 0; i < noStates; i++) {
                 for (int j = 0; j < noStates; j++){
                   denom =denom+ alfat[i][t]*transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1];
                 }
               }
               for (int i = 0; i < noStates; i++) {
                 gammat[i][t] = 0;
                 for (int j = 0; j < noStates; j++){
                   if (denom != 0){

                     diGammat[i][j][t] = (alfat[i][t]*transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1])/(denom+epsilon);
                     gammat[i][t]= gammat[i][t] +diGammat[i][j][t];
                   }
                 }
               }
             }
               //spec case for t =T-1
             denom = 0;
             for (int i = 0; i < noStates; i++) {
                 denom += alfat[i][observationSequence.length-1]; //*beta[observationSequence.length-1][i];
             }
             for (int i = 0; i < noStates; i++) {
               gammat[i][observationSequence.length-1] = (alfat[i][observationSequence.length-1])/(denom+epsilon); //*betat[i][observationSequence.length-1]
             }
           }

public static  void reEstimateInitialStateProb(){
    for (int i=0; i<noStates;i++){
        initialStateProb[i]=gammat[i][0];
    }
   }

   public static void reEstimateTransistionMatrix(int [] observationSequence){
       for (int i=0; i<noStates; i++){
           for (int j=0; j<noStates;j++){
               double numer=0;
               double denom =0;
               for (int t=0;t<observationSequence.length-1;t++){
                   numer=numer+diGammat[i][j][t];
                   denom=denom+gammat[i][t];
               }
               transitionMatrix[i][j]=numer/(denom+epsilon);
       }
     }
   }
        public static  void reEstimateEmissionMatrix(int [] observationSequence){
                for (int i=0; i< noStates;i++){
                        for (int j=0; j<noEmissions;j++){
                                double numer=0;
                                double denom=0;
                                for (int t=0; t<observationSequence.length;t++){
                                        if(observationSequence[t]==j){
                                                numer=numer+gammat[i][t];
                                        }
                                        denom=denom+gammat[i][t];
                                }
                                emissionMatrix[i][j]=numer/(denom+epsilon);

                        }
                }
        }

      public static  double computeConvergenceCondition2(int [] observationSequence){
          double logProb=0;
          for(int i=0;i<observationSequence.length;i++){
              logProb=logProb+Math.log(C[i]+epsilon);
          }

          logProb=-logProb;
          return logProb;
      }


  public static double BaumWelch(int [] observationSequence ){
    double[][] alfat = new double[noStates][observationSequence.length];
    double[][] betat = new double[noStates][observationSequence.length];
   diGammat = new double[noStates][noStates][observationSequence.length];
   gammat = new double[noStates][observationSequence.length];
   C= new double [observationSequence.length];

   denom =0;
   numer =0;

   double oldLogProb = -Double.MAX_VALUE;//-50000;//-Double.MAX_VALUE;//-100000000;
   int iters = 0;
   double logProb=0;
   int maxIters =100000000;  // 5
   boolean done = false;

   double diff = oldLogProb-logProb;
   //System.err.println("dif "+diff);
   alfat=alfapassScaling(observationSequence, alfat);
   betat =betapass(observationSequence,betat);
   diGammaGamma (betat,alfat,observationSequence);
   reEstimateInitialStateProb();
   reEstimateTransistionMatrix(observationSequence);
   reEstimateEmissionMatrix(observationSequence);
   logProb =computeConvergenceCondition2(observationSequence);

   iters += 1;



   while(logProb>oldLogProb && iters < maxIters  ){ //&& done==false){

     oldLogProb = logProb;
     //System.err.println("iters "+iters);

     alfat=alfapassScaling(observationSequence, alfat);
     betat =betapass(observationSequence,betat);
     diGammaGamma (betat,alfat,observationSequence);
     reEstimateInitialStateProb();
     reEstimateTransistionMatrix(observationSequence);
     reEstimateEmissionMatrix(observationSequence);
     logProb =computeConvergenceCondition2(observationSequence);

     iters += 1;
		 System.out.println("ITERATIONS: "+iters);


}
return logProb;
}


public static void main(String[] args){

	String fileName = args[0];
	HMM hmmobj = new HMM();
	File file = new File(fileName);

	Scanner input;

	try{
		input = new Scanner(file);
		String line4=input.nextLine();
		int [] observationSequence = observationFiller(line4);


		//System.out.println("OBSERVATION SEQUENCE"+Arrays.toString(observationSequence));
		double BM = BaumWelch(observationSequence);

		System.out.println("TRANSITION MATRIX"+Arrays.deepToString(hmmobj.transitionMatrix));
		System.out.println("EMISSION MATRIX"+Arrays.deepToString(hmmobj.emissionMatrix));
		System.out.println("INITIAL STATE MATRIX"+Arrays.toString(hmmobj.initialStateProb));

		}catch (FileNotFoundException exception) {
								System.out.println("The file does not exist");
							}

}
}
