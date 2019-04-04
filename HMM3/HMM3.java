import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.lang.Math;
import java.util.Collections;


public class HMM3{
  static double []c;

public static double [][] alfapass(double [][]transitionMatrix,double [][] emissionMatrix,int firstobsvalandindex, int noStates, int noTimeSteps, double []initialStateProb,int [] observationSequence){
  double []alfa0=new double [transitionMatrix[0].length];
  //double []
  c =new double[noTimeSteps];
  c[0]=0;

  for(int i=0;i<noStates;i++){
    alfa0[i]=initialStateProb[i]*emissionMatrix[i][firstobsvalandindex];
    c[0]=c[0]+alfa0[i];
  }
  //scale the alfa0(i)
  c[0]=1/c[0];
  for(int i=0;i<noStates;i++){
     alfa0[i]=alfa0[i]*c[0];
   }

  double [][]alfat=new double [noStates][noTimeSteps];

  for(int l=0; l<noStates;l++){
    alfat[l][0]=alfa0[l];

  }
  for (int t =1; t <noTimeSteps;t++){
    c[t]=0;
    for (int i=0;i<noStates;i++){
      for (int j=0; j<noStates;j++){
        alfat[i][t]=alfat[i][t]+alfat[j][t-1]*transitionMatrix[j][i];
      } // observationSequence[t] the value is the same as the index
      alfat[i][t]=alfat[i][t]*emissionMatrix[i][observationSequence[t]];
      c[t]=c[t]+alfat[i][t];
    }
    c[t]=1/c[t];
    for (int i=0; i<noStates;i++){
      alfat[i][t]=c[t]*alfat[i][t];
    }
}return alfat;
}
//
// public static double [][] betapass(double [][]transitionMatrix,double [][] emissionMatrix,int firstobsvalandindex, int noStates, int noTimeSteps, double []initialStateProb,int [] observationSequence){
//   double []beta0=new double [noStates];
//   double [][]betat=new double [noStates][noTimeSteps];
//   double tmp=0;
//   //System.out.println(Arrays.toString(c));
//   //System.out.println();
//   for(int i=0;i<noStates;i++){
//     tmp=c[noTimeSteps-1];
//     betat[i][noTimeSteps-1]=tmp;
//   }
//
//   for (int t =noTimeSteps-2; t <0;t--){
//     for (int i=0;i<noStates;i++){
//       for (int j=0; j<noStates;j++){
//         betat[i][t]=betat[i][t]+transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1];
//       } // observationSequence[t] the value is the same as the index
//       betat[i][t]=c[t]*betat[i][t];
//     }
// }return betat;
// }

public static double [][] betapass(double [][]transitionMatrix,double [][] emissionMatrix,int firstobsvalandindex, int noStates, int noTimeSteps, double []initialStateProb,int [] observationSequence){
  double []beta0=new double [noStates];
  double [][]betat=new double [noStates][noTimeSteps];
  double tmp=0;
  //System.out.println(Arrays.toString(c));
  //System.out.println();
  for(int i=0;i<noStates;i++){
    //tmp=c[noTimeSteps-1];
    betat[i][noTimeSteps-1]=c[noTimeSteps-1];
  }

  for (int t =noTimeSteps-2; t >=0;t--){
    for (int i=0;i<noStates;i++){
      for (int j=0; j<noStates;j++){
        betat[i][t]=betat[i][t]+transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1];
      } // observationSequence[t] the value is the same as the index
      betat[i][t]=c[t]*betat[i][t];
    }
}return betat;
}


public static double [][][] diGamma(double [][]transitionMatrix,double [][] emissionMatrix,int firstobsvalandindex, int noStates, int noTimeSteps, double []initialStateProb,int [] observationSequence, double [][]alfat,double [][]betat){

double [][][] diGammaMatrix=new double[noStates][noStates][noTimeSteps];
for (int t=0; t< noTimeSteps-1;t++){
  double denom = 0;
  for (int i=0; i< noStates;i++ ){
    for (int j=0; j< noStates;j++ ){
      denom=denom + alfat[i][t]*transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1];
  }
}
  for (int i=0; i<noStates; i++){
    for (int j=0; j<noStates; j++){
      diGammaMatrix[i][j][t]=(alfat[i][t]*transitionMatrix[i][j]*emissionMatrix[j][observationSequence[t+1]]*betat[j][t+1])/denom;
}
  }
}
return diGammaMatrix;
}

public static double [] reEstimateInitialStateProb(double[][]gamma,int noStates){
  double [] pi= new double [noStates];
  for (int i=0; i<noStates;i++){
    pi[i]=gamma[i][0];
  }
  return pi;
}

public static double [][] reEstimateTransistionMatrix(double[][]gammaMatrix,int noStates,double[][][]diGammaMatrix,int noTimeSteps){
  double [][] transitionMatrix=new double [noStates][noStates];
  for (int i=0; i<noStates; i++){
    for (int j=0; j<noStates;j++){
      double numer=0;
      double denom =0;
      for (int t=0;t<noTimeSteps-1;t++){
        numer=numer+diGammaMatrix[i][j][t];
        denom=denom+gammaMatrix[i][t];
      }
      transitionMatrix[i][j]=numer/denom;
    }
  }
return transitionMatrix;
}

public static double [][] reEstimateEmissionMatrix(double[][]gammaMatrix,int noStates,int noTimeSteps,int [] observationSequence,double [][] emissionMatrix){
  double [][]newEmissionMatrix=new double [noStates][emissionMatrix[0].length];
  for (int i=0; i< noStates;i++){
    for (int j=0; j<emissionMatrix[0].length;j++){
      double numer=0;
      double denom=0;
      for (int t=0; t<noTimeSteps;t++){
        if(observationSequence[t]==j){
          numer=numer+gammaMatrix[i][t];
        }
        denom=denom+gammaMatrix[i][t];
      }
      newEmissionMatrix[i][j]=numer/denom;
    }
  }
  return newEmissionMatrix;
}





public static double [][] gamma(double [][]transitionMatrix,double [][] emissionMatrix,int firstobsvalandindex, int noStates, int noTimeSteps, double []initialStateProb,int [] observationSequence, double [][]alfat,double [][]betat,double [][][] diGammaMatrix){

double [][] gammaMatrix=new double[noStates][noTimeSteps];
for (int t=0; t< noTimeSteps-1;t++){
  for (int i=0; i<noStates; i++){
    for (int j=0; j<noStates; j++){
      gammaMatrix[i][t]=gammaMatrix[i][t]+diGammaMatrix[i][j][t];
}
  }
}
double denom =0;
for (int i=0; i< noStates;i++){
  denom =denom + alfat[i][noTimeSteps-2];
}
for (int i=0; i< noStates;i++){
  gammaMatrix[i][noTimeSteps-2] =alfat[i][noTimeSteps-2]/denom;
}

return gammaMatrix;
}




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
  for(int i = 0; i < vectorF.length; i++) //Denna går igenom varje rad och för varje rad vill vi gå igen varje kolumn
  {
        vectorF[i]=Integer.parseInt(lineArray[itterator]);
        itterator++;

  }
  //String[] newArray = Arrays.copyOfRange(line, startIndex, endIndex);
  return vectorF;

}


public static double computeConvergenceCondition(int noTimeSteps){
  double logProb=0;
  for(int i=0;i<noTimeSteps;i++){
    logProb=logProb+Math.log(c[i]);
  }

  logProb=-logProb;
  return logProb;
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
  for(int r = 0; r < matrixF.length; r++) //Denna går igenom varje rad och för varje rad vill vi gå igen varje kolumn
  {
    //System.out.println("current row:"+r);
    for(int c = 0; c < matrixF[r].length; c++) // här går vi igenom varje element i en rad.
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

        String line=input.nextLine();
        double [][] transitionMatrix =matrixFiller(line);


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

        double [][] alfaf;
        double [][] betaf;
        double [][][] diGammat;
        double [][] gammat;
        //double [] pi;
        //System.out.println(transitionMatrix[0].length+" transmatrix");
        //System.out.println(emissionMatrix[0].length+"emissionMatrix");

        int maxIters=1000;
        int iters=0;
        double oldLogProb=Math.log(0);
        double logProb=0;
        //oldLogProb=logProb;
        alfaf=alfapass(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates,noTimeSteps,initialStateProb,observationSequence);
        betaf=betapass(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates,noTimeSteps,initialStateProb,observationSequence);
        diGammat=diGamma(transitionMatrix, emissionMatrix,firstobsvalandindex, noStates, noTimeSteps, initialStateProb, observationSequence,alfaf,betaf);

        //System.out.println(Arrays.deepToString(betaf));

        gammat=gamma(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates, noTimeSteps, initialStateProb,observationSequence,alfaf,betaf,diGammat);
        //System.out.println(Arrays.deegpToString(gammat));
        //System.out.println("betaf[0].length "+gammat[0].length);
        //System.out.println("alfaf.length "+gammat.length);

        initialStateProb=reEstimateInitialStateProb(gammat,noStates);
        //System.out.println(Arrays.toString(initialStateProb));
        transitionMatrix=reEstimateTransistionMatrix(gammat,noStates,diGammat,noTimeSteps);
        emissionMatrix=reEstimateEmissionMatrix(gammat,noStates,noTimeSteps,observationSequence,emissionMatrix);
        logProb=computeConvergenceCondition(noTimeSteps);

        iters++;

        while (iters<maxIters && logProb>oldLogProb){
          oldLogProb=logProb;
          alfaf=alfapass(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates,noTimeSteps,initialStateProb,observationSequence);
          betaf=betapass(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates,noTimeSteps,initialStateProb,observationSequence);
          diGammat=diGamma(transitionMatrix, emissionMatrix,firstobsvalandindex, noStates, noTimeSteps, initialStateProb, observationSequence,alfaf,betaf);
          gammat=gamma(transitionMatrix,emissionMatrix,firstobsvalandindex, noStates, noTimeSteps, initialStateProb,observationSequence,alfaf,betaf,diGammat);

          initialStateProb=reEstimateInitialStateProb(gammat,noStates);
          transitionMatrix=reEstimateTransistionMatrix(gammat,noStates,diGammat,noTimeSteps);
          emissionMatrix=reEstimateEmissionMatrix(gammat,noStates,noTimeSteps,observationSequence,emissionMatrix);
          logProb=computeConvergenceCondition(noTimeSteps);
          iters++;

        }
        //System.out.println("Iters: "+iters);


        //System.out.println(Arrays.deepToString(transitionMatrix)+"transmatrix");
        //System.out.println(Arrays.deepToString(emissionMatrix)+"emissionMatrix");


        String resTransM= Integer.toString(transitionMatrix.length)+" "+Integer.toString(transitionMatrix[0].length);
        System.out.print(resTransM);
        for(int i=0; i<transitionMatrix.length;i++){
          for (int j=0; j<transitionMatrix[0].length;j++){
            System.out.print(" ");
            System.out.print(transitionMatrix[i][j]);
          }

        }
        System.out.println("");
        String resEmMatrix= Integer.toString(emissionMatrix.length)+" "+Integer.toString(emissionMatrix[0].length);
        System.out.print(resEmMatrix);
        for(int i=0; i<emissionMatrix.length;i++){
          for (int j=0; j<emissionMatrix[0].length;j++){
            System.out.print(" ");
            System.out.print(emissionMatrix[i][j]);
          }

        }


                } catch (FileNotFoundException exception) {
                  System.out.println("The file does not exist");}
            }

        }
