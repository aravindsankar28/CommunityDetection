import java.util.ArrayList;
import java.util.HashMap;

/**
 * Contains helper functions for all mathematical computations such as likelihood, dot
 * product, matrix norms,etc.
 */
public class MathFunctions {

	static double heldOutLikelihood(ArrayList<ArrayList<Integer>> heldOutEdges, double[][] F, double[] W) {
		double likelihood = 0.0;
		int skippedEdges = 0;
		for (ArrayList<Integer> e : heldOutEdges) {
			int v1 = e.get(0);
			int v2 = e.get(1);
			double l = Math
					.log(1 - Math.exp(-MathFunctions.computePsi(v1, v2, F[v1], F[v2], W, Model.cosineSim, Model.alpha)));
			if (MathFunctions.dotProduct(F[v1], F[v2]) == 0 && Model.cosineSim[v1][v2] == 0)
				skippedEdges++;
			else
				likelihood += l;
		}
		return likelihood * (heldOutEdges.size() / (heldOutEdges.size() - skippedEdges));
	}

	static double likelihood(double[][] F, double[] W) {
		double Lg = 0.0;
		for (int u = 0; u < Model.G.size(); u++) {
			for (int v = u + 1; v < Model.G.size(); v++) {
				if (Model.G.get(u).containsKey(v))
					Lg += Math.log(
							1 - Math.exp(-MathFunctions.computePsi(u, v, F[u], F[v], W, Model.cosineSim, Model.alpha)));
				else
					Lg -= MathFunctions.computePsi(u, v, F[u], F[v], W, Model.cosineSim, Model.alpha);
			}
		}
		return Lg;
	}

	static double computePsi(int u, int v, double[] Fu, double[] Fv, double[] W, double[][] cosineSim, double alpha) {
		double val = 0.0;
		int C = W.length;
		for (int c = 0; c < C; c++) {
			val += Fu[c] * Fv[c] * W[c] * alpha;
		}
		val += (1.0 - alpha) * cosineSim[u][v];
		return val;
	}

	static double dotProduct(ArrayList<Double> x, ArrayList<Double> y) {
		double sum = 0.0;
		if (x.size() != y.size())
			return -1;

		for (int i = 0; i < x.size(); i++) {
			sum += x.get(i) * y.get(i);
		}
		return sum;
	}

	static double dotProduct(double[] x, double[] y) {
		double sum = 0.0;
		if (x.length != y.length) {
			System.out.println("error");
			return -1;
		}

		for (int i = 0; i < x.length; i++) {
			sum += x[i] * y[i];
		}
		return sum;
	}

	static double cosineSim(boolean[] x, boolean[] y) {
		int sum = 0;
		for (int i = 0; i < x.length; i++) {
			if (x[i] && y[i])
				sum++;
		}
		if (L2Norm(x) == 0 || L2Norm(y) == 0) {
			return 0.0;
		}
		return (sum / (L2Norm(x) * L2Norm(y)));
	}

	static double L2Norm(boolean[] x) {
		double sum = 0.0;
		for (boolean x1 : x) {
			if (x1)
				sum = sum + (1);
		}
		return Math.sqrt(sum);
	}

	static double L1Norm(HashMap<Integer, ArrayList<Double>> M) {
		double s = 0.0;
		for (int x : M.keySet()) {
			for (double y : M.get(x))
				s += Math.abs(y);
		}
		return s;
	}

	static double L1Norm(double[][] M) {
		double s = 0.0;
		for (int x = 0; x < M.length; x++) {
			for (int y = 0; y < M[x].length; y++)
				s += Math.abs(M[x][y]);
		}
		return s;
	}

	static double L2NormSq(HashMap<Integer, ArrayList<Double>> M) {
		double s = 0.0;
		for (int x : M.keySet()) {
			for (double y : M.get(x))
				s += y * y;
		}
		return s;
	}

	static double L2NormSq(double[][] M) {
		double s = 0.0;
		for (int x = 0; x < M.length; x++) {
			for (int y = 0; y < M[x].length; y++) {
				s += M[x][y] * M[x][y];
			}
		}
		return s;
	}
}
