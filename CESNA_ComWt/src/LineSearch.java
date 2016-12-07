
public class LineSearch {

	static double getLearningRateFU(int u, double[] updates, double[][] F, double[] W) {
		int C = updates.length;
		double beta = 0.5;
		double fx = likelihoodU(u, F[u], F, W);
		double t = 0.1;
		double controlParameter = 0.01;
		int iter = 0;
		boolean toContinue = true;
		do {
			iter++;
			if (iter > 10)
				return 0.0;
			double[] Fu_temp = new double[C];

			for (int c = 0; c < C; c++) {
				Fu_temp[c] = Math.max(0, F[u][c] + (t * updates[c]));
			}

			double fx_new = likelihoodU(u, Fu_temp, F, W);
			if (!Double.isFinite(fx_new)) {
				t = beta * t;
				continue;
			}
			if ((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * controlParameter) {
				t = beta * t;
				continue;
			} else
				toContinue = false;
		} while (toContinue);

		return t;
	}

	static double getLearningRateWC(double[] updates, double[][] F, double[] W) {
		double beta = 0.5;
		double fx = likelihood(F, W);
		double t = 0.1;
		double controlParameter = 0.01;
		int iter = 0;
		boolean toContinue = true;
		int C = updates.length;
		do {
			t = beta * t;
			iter++;
			if (iter > 10)
				return 0.0;
			double[] W_copy = new double[C];
			for (int u = 0; u < Model.G.size(); u++) {
				for (int c = 0; c < C; c++) {
					W_copy[c] = Math.max(0, W[c] + (t * updates[c]));
				}
			}
			double fx_new = likelihood(F, W_copy);
			if (((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * controlParameter)
					|| !Double.isFinite(fx_new)) {
				continue;
			} else
				break;

		} while (toContinue);

		return t;
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

	static double likelihoodU(int u, double[] Fu, double[][] F, double[] W) {
		double Lgu = 0.0;
		for (int v = 0; v < Model.G.size(); v++) {
			if (v == u)
				continue;
			if (Model.G.get(u).containsKey(v))
				Lgu += Math.log(1 - Math.exp(-MathFunctions.computePsi(u, v, Fu, F[v], W, Model.cosineSim, Model.alpha)));
			else
				Lgu -= MathFunctions.computePsi(u, v, Fu, F[v], W, Model.cosineSim, Model.alpha);
		}
		return Lgu;
	}

}
