import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.lang.Math;


public class HMM1{

public static double [][] matrixMultiplication(double [][] myMatrixA,double [][] myMatrixB ){

  int i, j, k, m, n;

  int rowsA = myMatrixA.length;
  int colsA = myMatrixA[0].length;
  int rowsB = myMatrixB.length;
  int colsB = myMatrixB[0].length;

  double[][] myMatrixC = new double [rowsA][colsB];


  //System.out.println("rader : "+rowsA);
  //System.out.println("kolumner "+colsB);

  for (i = 0; i < rowsA; i++) {

    for(j = 0; j < colsB; j++) {

     for(k = 0; k < colsA; k++) {
       //System.out.println(myMatrixA[i][k] * myMatrixB[k][j]);
      myMatrixC[i][j] += myMatrixA[i][k] * myMatrixB[k][j];

      //System.out.println(myMatrixA[i][k] * myMatrixB[k][j]);
     }
    }
  }return myMatrixC;

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



public static void main (String[] args){
        /**
        double [][] testmatrix1={{1,2,3},{2,2,2},{3,3,3}};
        double [][] testmatrix2={{1},{1},{1}};
        double [][] matrixres=matrixMultiplication(testmatrix1,testmatrix2);

        double[] cloneOfArray = Arrays.copyOf(matrixres[2], matrixres[2].length);
        System.out.println(Arrays.toString(cloneOfArray));
        **/



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

        //compute alfa0
        int firstobsvalandindex= observationSequence[0];
        //System.out.println(firstobsvalandindex);

        int noStates = transitionMatrix[0].length;
        int noTimeSteps =observationSequence.length;
        double c0 =0;
        double []alfa0=new double [transitionMatrix[0].length];
        for(int i=0;i<noStates;i++){
          alfa0[i]=initialStateProb[i]*emissionMatrix[i][firstobsvalandindex];
          c0=c0+alfa0[i];
          //System.out.println("emissionMatrix[i][firstobsvalandindex]: "+c0);
          //System.out.println("alfa0[i] "+alfa0[i]);
        }
        //System.out.println();


        //scale the alfa0(i)
        c0=1/c0;
        //tog bort detta skalning
        //for(int i=0;i<noStates;i++){
        //    alfa0[i]=alfa0[i]*c0;
            //System.out.println("alfa0[i] "+alfa0[i]);
        //  }

        //compute alfat(i)
        //double []alfa0=new double [transitionMatrix[0].length];
        double [][]alfat=new double [noStates][noTimeSteps];
        for(int l=0; l<alfa0.length;l++){
          alfat[l][0]=alfa0[l];

        }
        //System.out.println(Arrays.deepToString(alfat));
        Scanner reader = new Scanner(System.in);
        for (int t =1; t <noTimeSteps;t++){
          //System.out.println("New time step");
          double ct=0;
          for (int i=0;i<noStates;i++){
            //System.out.println("Current time step");
            for (int j=0; j<noStates;j++){
              //System.out.println("Previous time step");
              alfat[i][t]=alfat[i][t]+alfat[j][t-1]*transitionMatrix[j][i];
              //System.out.println("alfat[i][t]"+alfat[i][t]);
              //System.out.println("alfat[j][t-1]"+alfat[j][t-1]);
              //System.out.println("transitionMatrix[j][i]"+transitionMatrix[j][i]);
              //System.out.println("###########");
            } // observationSequence[t] the value is the same as the index
            alfat[i][t]=alfat[i][t]*emissionMatrix[i][observationSequence[t]];
            ct=ct+alfat[i][t];
          }
          //tog bort detta skalning
          //ct=1/ct;
          //for(int i=0;i<noStates;i++){
          //  alfat[i][t]=ct*alfat[i][t];
          //}

          //System.out.println();

        }
        double observationSequenceProbability=0;
        for( int i=0;i<noStates;i++){
          observationSequenceProbability =observationSequenceProbability+alfat[i][noTimeSteps-1];

        }
        System.out.println(observationSequenceProbability);









        } catch (FileNotFoundException exception) {
          System.out.println("The file does not exist");}
    }

}
