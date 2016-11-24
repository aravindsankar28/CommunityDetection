import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;

public class main {
	HashMap<Integer, ArrayList<Integer>> G;
	HashMap<Integer, ArrayList<Boolean>> X;
	HashMap<Integer, ArrayList<Double>> F; // u,c
	HashMap<Integer, ArrayList<Double>> W;// k,c
	int C; // No of communities
	int numAttr;
	double alpha;
	double eta;
	double lambda;
	int maxIter;
	double delta;
	double beta;

	public main(int com) {
		C = com;
		numAttr = 0;
		G = new HashMap<>();
		X = new HashMap<>();
		F = new HashMap<>();
		W = new HashMap<>();
		alpha = 0.5;
		eta = 0.001;
		lambda = 1;
		maxIter = 1000;
		beta = 0.7;
	}

	private void readGraph(String graphfilename, String attrfilename) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(graphfilename));
		String line;
		while ((line = br.readLine()) != null) {
			String[] edge = line.split(" ");
			int e0 = Integer.parseInt(edge[0]);
			int e1 = Integer.parseInt(edge[1]);
			if (!G.containsKey(e0)) {
				G.put(e0, new ArrayList<Integer>());
			}
			G.get(e0).add(e1);
		}
		br.close();
		br = new BufferedReader(new FileReader(attrfilename));
		while ((line = br.readLine()) != null) {
			String[] attributes = line.split(" ");
			if (numAttr == 0)
				numAttr = attributes.length - 1;
			int v = Integer.parseInt(attributes[0]);
			X.put(v, new ArrayList<Boolean>());
			for (int i = 1; i < attributes.length; i++) {
				X.get(v).add(Boolean.parseBoolean(attributes[i]));
			}
		}
		br.close();
		// adding back nodes that are in attrfile but not in graph
		for (int node : X.keySet()) {
			if (!G.containsKey(node)) {
				G.put(node, new ArrayList<Integer>());
			}
		}
		delta = -Math.log(1 - 1.0 / G.size());
	}

	private void initAffiliations() {
		// Compute the conductances of nodes
		HashMap<Integer, Double> conductance = new HashMap<>();
		for (int node : G.keySet()) {
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			neighbors.add(node);
			for (int n : G.get(node)) {
				neighbors.add(n);
			}
			int outedges = 0, inedges = 0;
			for (int v : neighbors) {
				for (int v1 : G.get(v)) {
					if (neighbors.contains(v1))
						inedges++;
					else
						outedges++;
				}
			}
			if (inedges == 0 && outedges == 0) {
				conductance.put(node, Double.MAX_VALUE);
			} else {
				conductance.put(node, (double) (outedges * 1.0 / inedges));
			}
		}

		conductance = (HashMap<Integer, Double>) mapUtil.sortByValue(conductance);
		ArrayList<Integer> alreadyAdded = new ArrayList<Integer>();
		int community = 0;
		double bias = 0.3;
		Random r = new Random();

		// Now initialize communities based on the conductance values computed
		for (int k : conductance.keySet()) {
			// if (conductance.get(k) < 0)
			// continue;
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k)) {
				neighbors.add(n);
			}
			for (int n : neighbors) {
				ArrayList<Double> Fu = new ArrayList<>();
				for (int c = 0; c < C; c++) {
					if (c == community) {
						Fu.add(0.7 + (0.3 * r.nextDouble()));
					} else {
						Fu.add(0.3 * r.nextDouble());
					}
				}
				F.put(n, Fu);
				alreadyAdded.add(n);
			}
			community++;
		}

		// Initialize weights W k,c
		for (int attr = 0; attr < numAttr; attr++) {
			W.put(attr, new ArrayList<Double>());
			for (int com = 0; com < C; com++) {
				W.get(attr).add((r.nextDouble() * 2.0) - 1.0);
			}
		}
	}

	// L_G
	double graphLikelihood(HashMap<Integer, ArrayList<Double>> F) {
		double Lg = 0.0;
		for (int u : G.keySet()) {
			for (int v : G.keySet()) {
				if (G.get(u).contains(v)) {
					Lg += Math.log(1 - Math.exp(-1 * MathFunctions.dotProduct(F.get(u), F.get(v))));
				} else
					Lg -= MathFunctions.dotProduct(F.get(u), F.get(v));
			}
		}

		return Lg;
	}

	// Lx
	double attributeLikelihood(HashMap<Integer, ArrayList<Double>> W) {
		double Lx = 0.0;
		for (int u : G.keySet()) {
			for (int k = 0; k < numAttr; k++) {
				boolean x = X.get(u).get(k);
				double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W.get(k), F.get(u))));
				int z = x ? 1 : 0;
				Lx += z * Math.log(Quk) + (1 - z) * Math.log(1 - Quk);
			}
		}
		return Lx;
	}

	double likelihood(HashMap<Integer, ArrayList<Double>> F, HashMap<Integer, ArrayList<Double>> W) {
		return (1 - alpha) * graphLikelihood(F) + alpha * attributeLikelihood(W) - lambda * MathFunctions.L1Norm(W);
	}

	double lineSearch(HashMap<Integer, ArrayList<Double>> updates) {
		double fx = (1 - alpha) * graphLikelihood(F) + (alpha) * attributeLikelihood(W);
		// System.out.println(fx);
		double t = 0.001;

		boolean toContinue = true;
		do {
			HashMap<Integer, ArrayList<Double>> F_copy = new HashMap<>(F);
			for (int u : G.keySet()) {
				for (int c = 0; c < C; c++) {
					F_copy.get(u).set(c, Math.max(0, F_copy.get(u).get(c) + (t * updates.get(u).get(c))));
				}
			}
			// System.out.println(t);
			// System.out.println(F_copy);
			double fx_new = (1 - alpha) * graphLikelihood(F_copy) + (alpha) * attributeLikelihood(W);
			//System.out.println("lg = "+graphLikelihood(F_copy));
		//	System.out.println("fx new = "+fx_new);
			// System.out.println(MathFunctions.L2NormSq(updates));
			System.out.println(fx_new-fx);
			
			if ((fx_new - fx) < t * MathFunctions.L2NormSq(updates) / 2.0) {
				t = beta * t;
			} else
				toContinue = false;
			

		} while (toContinue);
		return t;
	}

	void updateF() {
		// HashMap<Integer, ArrayList<Double>> newF = new HashMap<>();
		HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		double[] sumFvc = new double[C];

		for (int v : G.keySet()) {
			ArrayList<Double> Fv = F.get(v);
			for (int c = 0; c < C; c++) {
				sumFvc[c] += Fv.get(c);
			}
		}

		for (int u : G.keySet()) {
			// newF.put(u, new ArrayList<Double>());
			updates.put(u, new ArrayList<Double>());
			double[] derivative_G = new double[C];
			double[] derivative_X = new double[C];
			double[] diff = new double[C];
			for (int v : G.get(u)) {

				double dot = MathFunctions.dotProduct(F.get(u), F.get(v));
				double exp = Math.exp(-1 * dot);
				double fraction = exp / (1 - exp);

				for (int c = 0; c < C; c++) {
					derivative_G[c] += F.get(v).get(c) * fraction;
					diff[c] -= F.get(v).get(c);
					// derivative_G[c] -= F.get(v).get(c);
				}

			}
			for (int c = 0; c < C; c++) {
				diff[c] = diff[c] + sumFvc[c] - F.get(u).get(c);
				derivative_G[c] -= diff[c];
			}

			// derivative of Fu done

			for (int k = 0; k < numAttr; k++) {
				int z = X.get(u).get(k) ? 1 : 0;
				double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W.get(k), F.get(u))));
				double difference = (z - Quk);
				for (int c = 0; c < C; c++) {
					derivative_X[c] += difference * W.get(k).get(c);
				}
			}

			for (int c = 0; c < C; c++) {
				double update = (1 - alpha) * derivative_G[c] + (alpha) * derivative_X[c];
				updates.get(u).add(update);
			}
		}

		// eta = lineSearch(updates);
		// System.out.println(eta);
		for (int u : G.keySet()) {
			for (int c = 0; c < C; c++) {
				F.get(u).set(c,Math.max(0.0, F.get(u).get(c) + eta * updates.get(u).get(c)));
			}
		}
		// F = new HashMap<>(newF);
		// System.out.println(F);
	}

	void updateW() {
		HashMap<Integer, ArrayList<Double>> newW = new HashMap<>();
		for (int u : G.keySet()) {
			for (int k = 0; k < numAttr; k++) {
				double z = X.get(u).get(k) ? 1 : 0;
				double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W.get(k), F.get(u))));
				double diff = (z - Quk);
				newW.put(k, new ArrayList<Double>());
				for (int c = 0; c < C; c++) {
					double update = 0.0;
					update += diff * F.get(u).get(c);
					update *= alpha;
					if (W.get(k).get(c) > 0)
						update -= lambda;
					else
						update += lambda;
					newW.get(k).add(update * eta);
				}
			}
		}
	}

	void gradientAscent() {
		double startTime = System.currentTimeMillis();

		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = likelihood(F, W);
			System.out.println("Likelihood at iter " + iter + " " + likelihood);
			if (previousLikelihood != -1) {
				// System.out.println((likelihood - previousLikelihood) /
				// previousLikelihood);
				if ((previousLikelihood - likelihood) / previousLikelihood < 0.00001)
					break;
			}
			previousLikelihood = likelihood;

			updateF();
			updateW();
			System.out.println(eta);
			System.out.println("Time " + (System.currentTimeMillis() - startTime) / (1000.0 * (iter + 1)));
		}
	}

	void getCommunities() {
		ArrayList<ArrayList<Integer>> communities = new ArrayList<>();
		for (int c = 0; c < C; c++) {
			communities.add(new ArrayList<>());
		}

		for (int u : G.keySet()) {
			for (int c = 0; c < C; c++) {
				if (F.get(u).get(c) > delta) {
					communities.get(c).add(u);
				}
			}
		}

		for (ArrayList<Integer> comm : communities) {
			System.out.print("Circle");
			for (int x : comm)
				System.out.print(" " + x);
			System.out.println();
		}
	}

	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		main m = new main(10);
		String graphfilename = "facebook/0.edges";
		String attrfilename = "facebook/0.feat";
		m.readGraph(graphfilename, attrfilename);
		m.initAffiliations();
		System.out.println("Init done");
		m.gradientAscent();
		m.getCommunities();
	}

}
