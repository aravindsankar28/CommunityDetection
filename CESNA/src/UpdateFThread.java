class UpdateFThread extends Thread {

		// ArrayList<Integer> nodes;
		double[] sumFvc;
		// HashMap<Integer, ArrayList<Double>> updates;
		int startIndex;
		int endIndex;
		double[][] updates;
		double[][] F;
		double[][] W;

		UpdateFThread(double[] sumFvc, int startIndex, int endIndex, double[][] updates, double[][] F, double[][] W) {
			this.sumFvc = sumFvc;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
			this.updates = updates;
			this.F = F;
			this.W = W;
		}

		public void run() {
			updates = updateF(sumFvc, startIndex, endIndex, updates, F, W);
		}

		double[][] updateF(double[] sumFvc, int startIndex, int endIndex, double[][] updates, double[][] F, double[][] W) {
			int C = sumFvc.length;
			int numAttr = W.length;
			for (int u = startIndex; u <= endIndex; u++) {
				double[] derivative_G = new double[C];
				double[] derivative_X = new double[C];
				double[] diff = new double[C];
				for (int v : main.G.get(u).keySet()) {
					double dot = MathFunctions.dotProduct(F[u], F[v]);
					double exp = Math.exp(-1 * dot);
					double fraction = exp / (1 - exp);

					for (int c = 0; c < C; c++) {
						derivative_G[c] += F[v][c] * fraction;
						diff[c] -= F[v][c];

					}

				}
				for (int c = 0; c < C; c++) {
					diff[c] = diff[c] + sumFvc[c] - F[u][c];
					derivative_G[c] -= diff[c];
				}

				for (int k = 0; k < numAttr; k++) {
					int z = main.X[u][k] ? 1 : 0;
					double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
					double difference = (z - Quk);
					for (int c = 0; c < C; c++) {
						derivative_X[c] += difference * W[k][c];
					}
				}

				for (int c = 0; c < C; c++) {
					double update = (1 - main.alpha) * derivative_G[c] + (main.alpha) * derivative_X[c];
					updates[u][c] = update;

				}
			}

			return updates;
		}

	}