import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.lang.Math;


public class HMM0{

public static float [][] matrixMultiplication(float [][] myMatrixA,float [][] myMatrixB ){

  int i, j, k, m, n;

  int rowsA = myMatrixA.length;
  int colsA = myMatrixA[0].length;
  int rowsB = myMatrixB.length;
  int colsB = myMatrixB[0].length;

  float[][] myMatrixC = new float [rowsA][colsB];


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






public static float[][]  matrixFiller(String inputLine){
  String [] lineArray = inputLine.split("\\s+");  //text.split("\\s+");

  String[] cloneOfArray = Arrays.copyOf(lineArray, lineArray.length);
  //System.out.println(Arrays.toString(cloneOfArray));

  int rows,cols,itterator;
  rows=Integer.parseInt(lineArray[0]);
  cols=Integer.parseInt(lineArray[1]);
  //System.out.println(rows);
  //System.out.println(cols);

  itterator=2;

  float[][] matrixF=new float [rows][cols];

  //int i = startingNumber; i <= 100; i++
  for(int r = 0; r < matrixF.length; r++) //Denna går igenom varje rad och för varje rad vill vi gå igen varje kolumn
  {
    //System.out.println("current row:"+r);
    for(int c = 0; c < matrixF[r].length; c++) // här går vi igenom varje element i en rad.
    {
        //System.out.println("current col:"+c);
        //System.out.println(lineArray[itterator]);
        matrixF[r][c]=Float.parseFloat(lineArray[itterator]);
        itterator++;

    }
  }
  //String[] newArray = Arrays.copyOfRange(line, startIndex, endIndex);
  return matrixF;

}



public static void main (String[] args){
        /**
        float [][] testmatrix1={{1,2,3},{2,2,2},{3,3,3}};
        float [][] testmatrix2={{1},{1},{1}};
        float [][] matrixres=matrixMultiplication(testmatrix1,testmatrix2);

        float[] cloneOfArray = Arrays.copyOf(matrixres[2], matrixres[2].length);
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
        float [][] transitionMatrix =matrixFiller(line);

        //int[] cloneOfArray = Arrays.copyOf(transitionMatrix[2], transitionMatrix[2].length);
        //System.out.println(Arrays.toString(cloneOfArray));

        String line2=input.nextLine();
        float [][] emissionMatrix = matrixFiller(line2);

        String line3=input.nextLine();
        float [][] initialStateProb = matrixFiller(line3);

				System.out.println(Arrays.deepToString(initialStateProb));



        float[][]transitionMatrixInitalState=matrixMultiplication(initialStateProb, transitionMatrix);
        float [][]alfa2=matrixMultiplication(transitionMatrixInitalState, emissionMatrix);

        int final_row=alfa2.length;
        int final_col=alfa2[0].length;

        //System.out.println(Arrays.deepToString(alfa2));
        String res=Integer.toString(final_row)+" "+Integer.toString(final_col);
        //float tmp=100.00000001;
        for(float val: alfa2[0]){
          //System.out.println(val);
          //float tmp =Math.round(val);
          res+=" "+val;

        }
        System.out.println(res);


        } catch (FileNotFoundException exception) {
          System.out.println("The file does not exist");}
    }

}
