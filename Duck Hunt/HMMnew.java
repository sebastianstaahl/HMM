import java.util.*;
import java.util.Arrays;
import java.lang.Math;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.DoubleStream;
import static java.lang.System.*;
import java.util.Random;


public class HMMnew{
  static final Random noint = new Random();
  // Good seeds: 397169450 816905934 1033159171 1116287974 32

  //1483962880 1448742097  789487785 1019674898

//noint.nextInt()

  static final int rint=Math.abs(noint.nextInt());//2069414876);//1657773347);//1019674898      100p: 1657773347


  final double epsilon =(Double) Math.pow(10, -50);
  double []C;
  double [][] transitionMatrix;
  double [][] emissionMatrix;
  double [] initialStateProb;
  double [][][] diGammat;
  double [][] gammat;
  int noStates = 5;
  int noEmissions = 9;
  //static final int maxIters = 10000000; //TA BORT SEN
  double denom;
  double numer;

  public HMMnew() {
    //initilizeHmmModel();
    initilizeHmmModel3();
  }

// public void initilizeHmmModel2(){
//   this.transitionMatrix = new double[][] {{0.60, 0.10, 0.10, 0.10, 0.10}, {0.10, 0.60, 0.10, 0.10, 0.10}, {0.10, 0.10, 0.60, 0.10, 0.10}, {0.10, 0.10, 0.10, 0.60, 0.10}, {0.10, 0.10, 0.10, 0.10, 0.60}};
//   this.emissionMatrix = new double[][] {{0.050, 0.030, 0.050, 0.340, 0.200, 0.230, 0.030, 0.050, 0.040}, {0.400, 0.100, 0.450, 0.000, 0.050, 0.000, 0.000, 0.000, 0.000}, {0.000, 0.000, 0.000, 0.225, 0.050, 0.225, 0.050, 0.400, 0.050}, {0.120, 0.120, 0.120, 0.120, 0.040, 0.120, 0.120, 0.120, 0.120}, {0.125, 0.125, 0.125, 0.125, 0.000, 0.125, 0.125, 0.125, 0.125}};
//   this.initialStateProb = new double[] {0.21, 0.19, 0.17, 0.23, 0.20};
// }



public  double [][] alfaPassNO(int [] obs){
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

        public  double [][] matrixMultiplication(double [][] Matrix1,double [][] Matrix2 ){
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

public  double [] getColumn(double [][] matrix1,int c){
    double [] extractVector =new double[matrix1.length];
    for (int i=0; i<matrix1.length;i++){
        extractVector[i]=matrix1[i][c];

    }return extractVector;
}

public double alfaPassNOprob(int [] obs){
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



public void initilizeHmmModel(){
  this.transitionMatrix = new double[noStates][noStates];
  this.emissionMatrix = new double[noStates][noEmissions];
  this.initialStateProb = new double[noStates];

//.nextInt(upperbound-lowerbound) + lowerbound;
        //Random noint = new Random();
            Random no = new Random();
        //System.err.println("seed was: "+rint);


        //int rint=Math.abs(noint.nextInt());

            no.setSeed(rint); //3 bast
            Random r = new Random();
            r.setSeed(rint);
            //double randomValue = (8 + (13 - 8) * r.nextDouble());

            for (int i = 0; i < noStates; i++){
                double row_sum = 0;
                for (int j = 0; j < noStates; j++){
                    if(i==j){
                    this.transitionMatrix[i][j] = (no.nextDouble()+no.nextDouble()+no.nextDouble()+1)*40;
                    row_sum += this.transitionMatrix[i][j];
                }else{
                    this.transitionMatrix[i][j] = no.nextDouble()+1;
                    row_sum += this.transitionMatrix[i][j];
                }
                }
                for (int k = 0; k < noStates; k++){
                    this.transitionMatrix[i][k] = this.transitionMatrix[i][k]/(row_sum+epsilon);
                }
            }

            for (int i = 0; i < noStates; i++){
                double row_sum = 0;
                for (int j = 0; j < noEmissions; j++){
                    this.emissionMatrix[i][j] = no.nextDouble()+1;
                    row_sum += this.emissionMatrix[i][j];
                }
                for (int k = 0; k < noEmissions; k++){
                    this.emissionMatrix[i][k] = this.emissionMatrix[i][k]/(row_sum+epsilon);
                }
            }

            double row_sum = 0;
            for (int i = 0; i < noStates; i++){
                this.initialStateProb[i] = no.nextDouble()+1;
                row_sum += this.initialStateProb[i];
            }
            for (int k = 0; k < noStates; k++){
                this.initialStateProb[k] = this.initialStateProb[k]/(row_sum+epsilon);
            }

}



public void initilizeHmmModel3(){
  this.transitionMatrix = new double[noStates][noStates];
  this.emissionMatrix = new double[noStates][noEmissions];
  this.initialStateProb = new double[noStates];

//.nextInt(upperbound-lowerbound) + lowerbound;
        //Random noint = new Random();
            Random no = new Random();
        //System.err.println("seed was: "+rint);


        //int rint=Math.abs(noint.nextInt());

            //no.setSeed(rint); //3 bast  ///1516251283

            Random r = new Random();
            r.setSeed(1882442392);//rint);
            //double randomValue = (8 + (13 - 8) * r.nextDouble());

            for (int i = 0; i < noStates; i++){
                double row_sum = 0;
                for (int j = 0; j < noStates; j++){
                    if(i==j){
                    this.transitionMatrix[i][j] = 90;
                    row_sum += this.transitionMatrix[i][j];
                }else{
                    this.transitionMatrix[i][j] = (9 + (12 - 9)  * r.nextDouble());
                    row_sum += this.transitionMatrix[i][j];
                }
                }
                for (int k = 0; k < noStates; k++){
                    this.transitionMatrix[i][k] = this.transitionMatrix[i][k]/(row_sum+epsilon);
                }
            }

            for (int i = 0; i < noStates; i++){
                double row_sum = 0;
                for (int j = 0; j < noEmissions; j++){
                    this.emissionMatrix[i][j] = (9 + (12 - 9)  * r.nextDouble());
                    row_sum += this.emissionMatrix[i][j];
                }
                for (int k = 0; k < noEmissions; k++){
                    this.emissionMatrix[i][k] = this.emissionMatrix[i][k]/(row_sum+epsilon);
                }
            }

            double row_sum = 0;
            for (int i = 0; i < noStates; i++){
                this.initialStateProb[i] = (9 + (12 - 9) * r.nextDouble());
                row_sum += this.initialStateProb[i];
            }
            for (int k = 0; k < noStates; k++){
                this.initialStateProb[k] = this.initialStateProb[k]/(row_sum+epsilon);
            }


//System.err.println(Arrays.deepToString(this.transitionMatrix));
}


public  double [][] alfapassScaling(int [] observationSequence,double [][] alfat){

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


          public  double [][] betapass( int [] observationSequence ,double [][]betat){

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

           public void diGammaGamma (double [][] betat,double [][] alfat, int [] observationSequence){
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

public  void reEstimateInitialStateProb(){
    for (int i=0; i<noStates;i++){
        initialStateProb[i]=gammat[i][0];
    }
   }

   public void reEstimateTransistionMatrix(int [] observationSequence){
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
        public  void reEstimateEmissionMatrix(int [] observationSequence){
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

      public  double computeConvergenceCondition2(int [] observationSequence){
          double logProb=0;
          for(int i=0;i<observationSequence.length;i++){
              logProb=logProb+Math.log(C[i]+epsilon);
          }

          logProb=-logProb;
          return logProb;
      }


  public double BaumWelch(int [] observationSequence ){
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
   int maxIters =6;  // 5
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


}
return logProb;
}
}
