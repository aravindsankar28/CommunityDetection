import java.util.ArrayList;

class UpdateWThread extends Thread {

		ArrayList<Integer> attributes;
		
		double[][] updates;
		int startIndex;
		int endIndex;
		double[][]F;
		double[][] W;

		UpdateWThread(int startIndex, int endIndex, double[][] updates, double[][]F, double[][] W) {
			this.startIndex = startIndex;
			this.endIndex = endIndex;
			this.updates = updates;
			this.F = F;
			this.W = W;
		}


		double[][] updateW(int startIndex, int endIndex, double[][] updates, double[][]F, double[][] W) {
			int C = W[0].length;
			for (int k = startIndex; k <= endIndex; k++) {
				for (int c = 0; c < C; c++) {
					double update = 0.0;
					for (int u = 0; u < main.G.size(); u++) {
						double z = main.X[u][k] ? 1 : 0;
						double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
						double diff = (z - Quk);
						update += diff * F[u][c];
					}
					update *= main.alpha;
					if (W[k][c] > 0)
						update -= main.lambda;
					else
						update += main.lambda;
					updates[k][c] = update;
				}
			}
			return updates;
		}
		
		public void run() {
			updates = updateW(startIndex, endIndex, updates, F, W);
		}
	}

