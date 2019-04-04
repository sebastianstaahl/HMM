//kolla hur manga modeller vi sparar. om vi bara sparar 6 stycken alltid.
//lagg till en vektor med logprobs for varje modell och jamfor om det finns ngn battre

/*
1. valj den basta modellen for varje fagel ex pigeon no 20. Gor detta for alla faglar
2.valj fageln med bast P(O|lambda)

*/

import java.util.*;
import java.util.stream.*;
import java.util.Arrays;
import java.util.Arrays.*;
import java.util.Iterator;

class Player {
  int time = 0;
  //public HMMnew [] birdHmmList=new HMMnew[Constants.COUNT_SPECIES];
  public ArrayList[] birdModelsHMM=new ArrayList[6]; //We wanna save all HMM models for every bird. Therefore we create a dynamic matrix
  final double epsilon =(Double) Math.pow(10, -50);
  int noStates = 5;
  int noEmissions = 9;
  int cG, fG,cS,tS;
	double guessLimit = -1000;
	double shootLimit = 460;


  public HMMnew [] hmmsPerSpec =new HMMnew[Constants.COUNT_SPECIES];

  int [] lGuess;


    public Player() {

      cG = 0;
        fG = 0;
      cS=0;
      tS=0;

      for (int j=0; j<Constants.COUNT_SPECIES;j++){

        //birdModelsHMM[j]= new ArrayList(); //we wanna add a dynamic arraylist at
        birdModelsHMM[j]= new ArrayList<HMMnew>();
        //every position in arraylist. but we want to outer list to be an array
        //so that we can use the itterator
      }

    }


     public  double [] getColumn(double [][] matrix1,int c){
                double [] extractVector =new double[matrix1.length];
                for (int i=0; i<matrix1.length;i++){
                        extractVector[i]=matrix1[i][c];

                }return extractVector;
        }

		public double [] findShootEmission(double [] emissionProbs,int [] pMovements){
		  double [] res =new double [noEmissions];
		  int [] noOcurr =new int [noEmissions];


		  for(int i=0;i<pMovements.length;i++){

		    if (pMovements[i]!=-1){ //ev ta bort detta
		    noOcurr[pMovements[i]]+=1;
		    res[pMovements[i]]+=emissionProbs[i];
		  }else{
		    res[9]+=emissionProbs[i];
		    noOcurr[9]+=1;
		  }
		  }
		  for(int i=0;i<res.length;i++){
		    if(res[i]!=0){
		    res[i]=res[i]/noOcurr[i];
		  }
		}
		return res;

		}

		public int findBestBirdToShoot(double [] logProbsBirds ,int [] bestBirdsSpieces){
		  int index = indexOfMaxFinderlog(logProbsBirds);

		  for(int i =0; i <bestBirdsSpieces.length ; i++){
		    if (bestBirdsSpieces[index] == Constants.SPECIES_BLACK_STORK || bestBirdsSpieces[index] == Constants.SPECIES_UNKNOWN){
		    logProbsBirds[index]=0;}
		  }

		  index = indexOfMaxFinderlog(logProbsBirds);

		  return index;

		}

		/**
		 * Shoot!
		 *
		 * This is the function where you start your work.
		 *
		 * You will receive a variable pState, which contains information about all
		 * birds, both dead and alive. Each bird contains all past moves.
		 *
		 * The state also contains the scores for all players and the number of
		 * time steps elapsed since the last time this function was called.
		 *
		 * @param pState the GameState object with observations etc
		 * @param pDue time before which we must have returned
		 * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
		 */
		//################################
		public Action shoot(GameState pState, Deadline pDue) {
			time++;
			//Wait until the 88 timestep until we shoot.
			int villkor =86+99*pState.getRound();
		  if( pState.getRound()>2 && time>villkor){ // best 2
				//Inital max value of logprob.
		    double maxLogProb = -Double.MAX_VALUE;
				//Create 3 lists
				//one with the probability for the bird specie which was best, for each bird in the environment.
				//One for the hmm models, one hmm model for each of the birds in the environment.
				//and one that contains which specie had the highest probability for each bird in the environment.
		    double [] logProbsBirds = new double [pState.getNumBirds()];
		    HMMnew [] bestBirdsHMMs = new HMMnew [pState.getNumBirds()];
		    int [] bestBirdsSpieces = new int [pState.getNumBirds()];
				//Loop through the birds in the environment
		    for(int i = 0; i < pState.getNumBirds(); ++i){
		      if(pState.getBird(i).isAlive() && pState.getBird(i).getSeqLength()>63){
						int noBirdObs = pState.getBird(i).getSeqLength();
			      int [] observstionSeq = new int [noBirdObs];
			      for (int p=0; p < noBirdObs; p++){
			        observstionSeq[p]=pState.getBird(i).getObservation(p);
			      }
						//Loop through all of the species
			      for(int g = 0; g<Constants.COUNT_SPECIES; g++){
							Iterator<HMMnew> it;
			      	it = birdModelsHMM[g].iterator();
			      	//for (HMMnew birdHmm: birdHmmList){
			      	ArrayList<Double> lopProbBirdModelPicker =new ArrayList <Double>();
							//Loop through all of the models for the given specie.
				      while (it.hasNext()){
				        HMMnew birdHmm =it.next();
								//Retrieve the logprob from alfapass.
				        double probObsSeqLog =birdHmm.alfaPassNOprob(observstionSeq);//alfapassnoscaling
								//if the model generates a higher maxprobability, save relevant information.
								//We will handle the cases where Stork is the best below (row 162).
				        if(probObsSeqLog>maxLogProb){
				          maxLogProb=probObsSeqLog;
				          logProbsBirds[i]=probObsSeqLog;
				          bestBirdsSpieces[i] = g;
				          bestBirdsHMMs[i]=birdHmm;

		              }
		            }
		        }
		      }
		    }
				//Check if the bird with the highest probability is a stork
				//If it is a stork return the second best bird to shoot, out of all of the
				//birds in the environment.
				int bestBird = findBestBirdToShoot(logProbsBirds, bestBirdsSpieces);
				//Get the specie for the best bird.
		    int theBestBirdsSpieces=bestBirdsSpieces[bestBird];
				//Number of observations for the best bird.
		    int noBirdObs = pState.getBird(bestBird).getSeqLength();
				//A list for the observations for the bist bird.
		    int [] observstionSeq = new int [noBirdObs];
				//If the probability is higher than a criteria and best bird to shoot isn't unknown specie or stork(double check)
		    if(maxLogProb> shootLimit && bestBird !=-1 && theBestBirdsSpieces != Constants.SPECIES_BLACK_STORK){
					//Create an iterator for the best specie.
		      Iterator<HMMnew> it2;
					it2 = birdModelsHMM[theBestBirdsSpieces].iterator();

					//Fill the observation list for the best bird.
		      for (int p=0; p < noBirdObs; p++){
		        observstionSeq[p]=pState.getBird(bestBird).getObservation(p);
		      }

		      //Here we are going to loop through all of the HMM-models for the best specie.
					//maxnextemission is the index to state with the highest next emission probability.
		      int maxnextemission=0;
					//maxprobnextemission is a variable that will containt the highest probability of one state.
		      double maxprobnextemission= -Double.MAX_VALUE;

					//Loop through the hmm-models.
		      while(it2.hasNext()){
		      HMMnew birdHmm = it2.next();//bestBirdsHMMs[bestBird];//
					//Calculate logprob.
		      double [][] birdAlfaPass =birdHmm.alfaPassNO(observstionSeq);
					//Get the last column.
		      double [] StateProbLastTimeStep =birdHmm.getColumn(birdAlfaPass,birdAlfaPass[0].length-1);
					//Need to convert it to a matrix for our matrixMultiplication
		      double [][] StateProbLastTimeStepMatrix = new double [1][StateProbLastTimeStep.length];
					//Insert the values into the matrix.
		      for(int d=0;d<StateProbLastTimeStepMatrix.length;d++){
		        for(int k=0;k<StateProbLastTimeStepMatrix[d].length;k++){
		          StateProbLastTimeStepMatrix[d][k]=StateProbLastTimeStep[k];
		      }
		    }
					//Multiplication of State probability distribution x TransitionMatrix
		      double [][] probNextState = birdHmm.matrixMultiplication(StateProbLastTimeStepMatrix,birdHmm.transitionMatrix);
					//Multiplying the results with EmissionMatrix to get next emission distribution.
		      double [][] emissionDistributionNextState = birdHmm.matrixMultiplication(probNextState ,birdHmm.emissionMatrix);
					//Get the vector from the matrix.
		      double [] emissionDistributionNextStateVector=emissionDistributionNextState[0];
					//Find the most probable next emission.
		      double probnextemission=maxFinder(emissionDistributionNextStateVector);
					//Find the index to the most probable next emission.
		      int nextemission = indexOfMaxFinder(emissionDistributionNextStateVector);
					//Check if it is the biggest probability.
			    if(probnextemission>maxprobnextemission){
			      maxprobnextemission=probnextemission;
			      maxnextemission=nextemission;
			     }
				 }
				 //Chech if it fulfulls the criteria
		      if (maxprobnextemission>0.61){ //&& sumBirdEmissions[Constants.SPECIES_BLACK_STORK]<3 ){//probMax>13 && emind != Constants.SPECIES_BLACK_STORK && sumBirdEmissions[Constants.SPECIES_BLACK_STORK]<3){ //lagg till vilkor typ: om blackstork prob ar under
		        tS++;
						//Return the action.
						//Bestbird is which of the birds in the environment we want to shoot.
						//Maxnextemission is in which direction we are predicting the bird is going to fly.
		        return new Action(bestBird,maxnextemission);
		      	}else{
		        	return cDontShoot;
		      	}
	}
		// for the first round or default in other rounds when u don't have time left or the prob for an emission was to low
		}
		 return cDontShoot;
		}



    public  double maxFinder(double [] vector1){
        double biggest=vector1[0];
        for(int i=1; i<vector1.length;i++){
            if(vector1[i]>biggest){
                biggest=vector1[i];
            }
        } return biggest;
    }

    public  int indexOfMaxFinder(double [] vector1){
                double biggest=vector1[0];
                int indexOfBiggest=0;
                for(int i=1; i<vector1.length;i++){
                        if(vector1[i]>biggest){
                                biggest=vector1[i];
                                indexOfBiggest=i;
                        }
                } return indexOfBiggest;
        }

    public  int indexOfMaxArrayFinder(ArrayList<Double> AL){
                double biggest=AL.get(0);
                int indexOfBiggest=0;

                for(int i=1; i<AL.size();i++){
                        if(AL.get(i)>biggest){
                                biggest=AL.get(i);
                                indexOfBiggest=i;
                        }
                } return indexOfBiggest;
        }

    public  double maxArrayListFinder(ArrayList<Double> AL){
        double biggest=AL.get(0);
        for(int i=1; i<AL.size();i++){
            if(AL.get(i)>biggest){
                biggest=AL.get(i);
            }
        } return biggest;
    }
    public  double maxFinderlog(double [] vector1){
        double biggest=-10000000;//vector1[0];

        for(int i=0; i<vector1.length;i++){
            if(vector1[i]>biggest && vector1[i]!=0){
                biggest=vector1[i];
            }
        } return biggest;
    }

    public  int indexOfMaxFinderlog(double [] vector1){
        double biggest=-Double.MAX_VALUE;//vector1[0];
        int indexOfBiggest=-1;

        for(int i=0; i<vector1.length;i++){
            if(vector1[i]>biggest && vector1[i]!=0){
                biggest=vector1[i];
                indexOfBiggest=i;
            }
        }if(biggest==-Double.MAX_VALUE){
          return -1;
        }
    return indexOfBiggest;
      }


     public double secondLargestValue(double[] vector ) {
            double l = vector[0];
            double sl = vector[0];
            for (int i = 0; i < vector.length; i++) {
                if (vector[i] > l) {
                    sl = l;
                    l = vector[i];
                } else if (vector[i] > sl) {
                    sl= vector[i];
                }
            }
    return sl;
    }

    public double secondLargestValuelog(double[] vector ) {
       double l = -10000000;//vector[0];
       double sl = -100000000;//vector[0];
       for (int i = 0; i < vector.length; i++) {
         if (vector[i] > l && vector[i]!=0) {
           sl = l;
           l = vector[i];
         } else if (vector[i] > sl && vector[i]!=0) {
           sl= vector[i];
         }
       }
   return sl;
   }



    public int probabilityComparer(double [] modelSelectionprobabilities){
      //System.err.println(Arrays.toString(modelSelectionprobabilities));

      double maxProb = maxFinderlog(modelSelectionprobabilities);
      int maxIndex = indexOfMaxFinderlog(modelSelectionprobabilities);

      double secondLargestProb =  secondLargestValuelog(modelSelectionprobabilities);
      int secLarIndex = Arrays.asList(modelSelectionprobabilities).indexOf(secondLargestProb);

      double div =maxIndex/(secLarIndex+epsilon);
      double minProb =-420;//(Double) Math.pow(10, -60);

      //Check taht they are not to close
      if(maxProb>minProb){ //om vi kan fa tbx vardena i riktiga sannolikehter da vill vi aven kolla att sannolikehen ar stor nog
        return maxIndex;
      }else{
        return Constants.SPECIES_UNKNOWN;

      }
    }

		/**
		 * Guess the species!
		 * This function will be called at the end of each round, to give you
		 * a chance to identify the species of the birds for extra points.
		 *
		 * Fill the vector with guesses for the all birds.
		 * Use SPECIES_UNKNOWN to avoid guessing.
		 *
		 * @param pState the GameState object with observations etc
		 * @param pDue time before which we must have returned
		 * @return a vector with guesses for all the birds
		 */
    public int[] guess(GameState pState, Deadline pDue) {
      /*
       * Here you should write your clever algorithms to guess the species of
       * each bird. This skeleton makes no guesses, better safe than sorry!
       */
      Random rand = new Random();
			//Our guess array.
      lGuess = new int[pState.getNumBirds()];

			//Fill the array with specie unknown.
      for (int i = 0; i < pState.getNumBirds(); ++i)
          lGuess[i] = Constants.SPECIES_UNKNOWN; //ger ingen gissning

			//If we don't observe any birds we dont want to make any guesses.
      if (pState.getNumBirds()==0){
            return lGuess;
          }else{
						//If it is the first round, we guess on random birds.
						//we want to guess to even though we dont know anything in round 0,
						//because we need to get HMM-models.
            if(pState.getRound()==0){
              for(int i = 0; i < pState.getNumBirds(); ++i){
                  lGuess[i] = rand.nextInt(Constants.COUNT_SPECIES-1); //change to the bird we wanna guess for each index
                }

								//If it's not the first round we want to use our HMM-models to make the guesses.
           }else{
						 //We loop through all of the birds in the environment.
             for(int j = 0; j < pState.getNumBirds(); ++j){
							 //Create a list for each bird to save the observation sequence.
               int noBirdObs = pState.getBird(j).getSeqLength();
               int [] observstionSeq = new int [noBirdObs];
							 //Save the observations in the list.
               for (int p=0; p < noBirdObs; p++){
                 //Check that the bird is alive before saving the observation.
                 if(pState.getBird(j).getObservation(p)!=Constants.MOVE_DEAD){
                 observstionSeq[p]=pState.getBird(j).getObservation(p);
               }
               }
							 //Initail value for the convergence check.
               double maxProb =-100000;
							 //index to the specie that gives the highest probability.
               int maxIndex=0;
							 //Iterate through the different species.
               for(int d = 0; d<Constants.COUNT_SPECIES; d++){

								//Create an iterator for each specie.
               Iterator<HMMnew> iterator;
               iterator = birdModelsHMM[d].iterator();
							 //Iterate through all of the HMM-models for the specie.
                while (iterator.hasNext()){ //we wanna look at all the models. which we save after every round.
                  HMMnew hmmBirdSpiecesModelround =iterator.next();
									//Call on alfapass.
                  double probObsSeq = hmmBirdSpiecesModelround.alfaPassNOprob(observstionSeq);
									//Check if the logprog value is greater than the currently max logprob.
                  if (probObsSeq>maxProb){
                    maxProb=probObsSeq;
                    maxIndex=d;
                  }
              }
            }	//Check if the maxprob is greater than a critera before we guess.
						//If we are in the three first rounds, we wanna guess if the maxprob isn't greater because we want more HMM-models.
              if(maxProb > guessLimit || pState.getRound() < 3){
              lGuess[j]=maxIndex;//=maxin;//indexOfMaxFinder(meanProbperSpieces);//maxIndex;//indexOfMaxFinder(meanProbperSpieces);//maxIndex; //birdGuess;
            }

           }

         }
          }
      return lGuess;
  }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
       cS++;
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
      int noBirds=pSpecies.length;
			//
      for(int i=0; i<noBirds;i++){
				//for each bird we want to create an HMM-model that we will train with
				//the observation sequence that we receive

				//Check that the specie is know so we dont train a modell for an unknown
				//specie.
        if(pSpecies[i] != Constants.SPECIES_UNKNOWN){

				//Create a list for the observation sequence.
        int numOfObs = pState.getBird(i).getSeqLength();
        int[] observationSequence = new int[numOfObs];
				//Fill the list with the observations of the bird.
        for (int k = 0; k < numOfObs; k++) {
          if(pState.getBird(i).getObservation(k)!=Constants.MOVE_DEAD){
            observationSequence[k] = pState.getBird(i).getObservation(k);
          }
        }

				//Create an HMM-model for the bird. (one for each bird).
        HMMnew birdHMM = new HMMnew();
				//Train the HMM-modell using Baum Welch.
        double logProb=birdHMM.BaumWelch(observationSequence);
				//Add the HMM-model to the same specie as the bird currently observing.
				//Adding to an ArrayList where each index is the bird specie, and at
				//each index there is an arraylist containing the HMM-models for the specie.
        birdModelsHMM[pSpecies[i]].add(birdHMM);



    if(lGuess[i] == pSpecies[i]){
      cG=cG+1;
    } else {
      fG=fG+1;
    }
  }

    int diff =tS-cS;


  }

}
    public static final Action cDontShoot = new Action(-1, -1);
}
