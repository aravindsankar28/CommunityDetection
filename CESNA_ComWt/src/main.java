import java.util.HashMap;

/**
 * This program calls our model multiple times to determine
 */
public class Main {
	static Model model;
	static int nRandLim = 3;

	// By running our model multiple times, this function returns the best
	// number of communities.
	static int chooseCommunites(String graphFileName, String attrFileName, String outputFilename) throws Exception {
		HashMap<Integer, Double> perf = new HashMap<>();
		int bestK = 0;
		double bestKLikelihood = (double) 0;
		int lowerLimit = 4;
		int upperLimit = 8;
		if (model.V >= 400) {
			lowerLimit = 6;
			upperLimit = 14;
		}
		for (int K = lowerLimit; K <= upperLimit; K += 2) {
			System.out.println(K);
			model.C = K;
			double avgL = 0.0;
			for (int nRand = 0; nRand < nRandLim; nRand++) {
				System.out.println(nRand);
				model.readGraph(graphFileName, attrFileName);
				model.performHoldOut(true);
				model.driver(graphFileName, attrFileName, outputFilename, false);
				double L = MathFunctions.heldOutLikelihood(model.heldOutEdges, model.F, model.W);
				avgL += L;
				model.performHoldOut(false);
			}
			avgL /= nRandLim;
			perf.put(K, avgL);
			if (bestKLikelihood == 0.0) {
				bestKLikelihood = avgL;
				bestK = K;
			} else {
				if (avgL > bestKLikelihood) {
					bestKLikelihood = avgL;
					bestK = K;
				}
			}
		}
		return bestK;
	}

	public static void main(String[] args) throws Exception {
		model = new Model(0);
		int graphNum = 3980;

		String graphFileName = "facebook/" + String.valueOf(graphNum) + ".edges";
		String attrFileName = "facebook/" + String.valueOf(graphNum) + ".feat";
		String outputFilename = "circles.txt";
		// Change to command line args.
		
		model.readGraph(graphFileName, attrFileName);
		int bestK = chooseCommunites(graphFileName, attrFileName, outputFilename);
		model.C = bestK;
		model.stoppingCondition /= 1000;
		model.readGraph(graphFileName, attrFileName);
		model.driver(graphFileName, attrFileName, outputFilename,  true);
		System.out.println("-");
		System.out.println(bestK);
	}
}
