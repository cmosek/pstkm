package swnd;

import ilog.concert.*;
import ilog.cplex.*;

public class Model {
    //Wersja V4

    int d_len = 6;
    int v_len = 4;
    int p_len = 12;

    //Wersja V3
    /*
    int d_len = 6;
    int v_len = 3;
    int p_len = 12;

    */
        //Wersja V2
    /*
        int v_len = 2;
        int d_len = 2;
        int p_len = 2;
    */

    int version = 4;

    //	Zbiory
    Integer[] D = new Integer[d_len];    //	Zbior zapotrzebowan
    Integer[] V = new Integer[v_len];    //	Zbior wezlow
    Path[] P = new Path[p_len];        //	Zbior sciezek


    // Wartosci stale
    Integer[][][] delta_vdp = new Integer[V.length][D.length][P.length];    // [v][d][p], {0,1}	=1-jezeli wezel v lezy na sciezce p zapotrzebowania d,=0-w innym przyp
    Integer[][] s_vd = new Integer[V.length][D.length];            // [v][d], {0,1} =0 –wezel v jest wezlem startowym dla zapotrzebowania d,=1-w innym przyp
    Integer[][] t_vd = new Integer[V.length][D.length];            // [v][d], {0,1} =0 –wezel v jest wezlem koncowym dla zapotrzebowania d,=1-w innym przyp
    Integer[] T_v = new Integer[V.length];                // [v], maksymalna szybkosc z jaka moze wysylac dane wezel, T_v >= 0
    Integer[] R_v = new Integer[V.length];                // [v], maksymalna szybkosc z jaka moze odbierac dane wezel, R_v >= 0
    Double[] pc_v = new Double[V.length];            // [v], wsp. zuzycia energii wezla v
    ///TODO zdefiniowac r


    public Model(double r) {
        model(r);
    }

    public void model(double r) {
        try {
            IloCplex cplex = new IloCplex();


            startData();
            fill_delta_s_t();


            // Zmienne
            /// Definicja zmiennej x_dp, ruch zapotrzebowania d na sciezce p, x_dp >= 0
            IloIntVar[][] x_dp = new IloIntVar[D.length][P.length];
            for (Integer d : D) {
                for (Path p : P) {

                    x_dp[d - 1][p.getNr() - 1] = cplex.intVar(0, Integer.MAX_VALUE, "x_" + d + "-" + p.getNr());
                }
            }

            /// Definicja zmiennej h_d, volumen zapotrzebowania d, h_d >= 0
            IloIntVar h_d = cplex.intVar(0, Integer.MAX_VALUE, "h_d");


            //System.out.println(delta_vdp[0]);
            for (Integer v : V) {
                for (Integer d : D) {

                    int[] ret = new int[p_len];
                    for (int i = 0; i < ret.length; i++)
                        ret[i] = delta_vdp[v - 1][d - 1][i];

                    cplex.addLe(cplex.sumx_dp[d - 1], h_d);
                    //cplex.addEq(cplex.scalProd(x_dp[d - 1], ret), h_d);
                }
            }

			/*
            // Ograniczenia

			IloIntExpr sum;

			//for (Integer d : D) {
                double[] objvals = {1.0, 0};
                //cplex.addMaximize();
				//sum = cplex.sum(x_dp[d-1]);
				cplex.addEq(cplex.scalProd(x_dp[0], objvals), h_d);   // Moze trzeba bedzie to inaczej zapisac
			//}
            double[] objvals2 = {0, 1.0};
            cplex.addEq(cplex.scalProd(x_dp[1], objvals2), h_d);
			*/

            //TODO wywalić
            //cplex.addGe(h_d, 0.0);
            //cplex.addLe(h_d, Integer.MAX_VALUE);


            // Ograniczenie na maksymalną szybkość wysyłania danych przez węzeł
            for (Integer v : V) {
                IloLinearIntExpr lhs = cplex.linearIntExpr();
                for (Integer d : D) {
                    for (Path p : P) {
                        Integer tmp = t_vd[v - 1][d - 1];
                        Integer temp = delta_vdp[v - 1][d - 1][p.getNr() - 1];
                        IloIntVar x = x_dp[d - 1][p.getNr() - 1];
                        //System.out.print(v.toString());
                        //System.out.println(x);
                        lhs.addTerm(temp * tmp, x);
                    }
                }
                cplex.addLe(lhs, T_v[v - 1] * r);
            }


            /// Ograniczenie na maksymalną szybkość odbierania danych przez węzeł

            for (Integer v : V) {
                IloLinearIntExpr lhs = cplex.linearIntExpr();
                for (Integer d : D) {
                    for (Path p : P) {
                        Integer temp = delta_vdp[v - 1][d - 1][p.getNr() - 1] * s_vd[v - 1][d - 1];
                        lhs.addTerm(temp, x_dp[d - 1][p.getNr() - 1]);
                    }
                }
                cplex.addLe(lhs, R_v[v - 1]);
            }


            /// Ograniczenie na maksymalne zużycie energii przez węzeł

            for (Integer v : V) {
                IloLinearIntExpr lhs = cplex.linearIntExpr();
                for (Integer d : D) {
                    for (Path p : P) {
                        Integer temp = delta_vdp[v - 1][d - 1][p.getNr() - 1] * (s_vd[v - 1][d - 1] + t_vd[v - 1][d - 1]);
                        lhs.addTerm(temp, x_dp[d - 1][p.getNr() - 1]);
                    }
                }
                Double d = ((double) T_v[v - 1] * r + (double) R_v[v - 1]) * pc_v[v - 1];
                Integer rhs = d.intValue();
                cplex.addLe(lhs, rhs);
            }


            // Cel
            cplex.addMaximize(h_d);
// turn off the presolver
            //cplex.setParam(IloCplex.BooleanParam.PreInd, false);
// select the primal simplex method
            //cplex.setParam(IloCplex.IntParam.RootAlg,IloCplex.Algorithm.Primal);
            if (cplex.solve()) {
                //String s = "asd.lp";
                //cplex.exportModel(s);
                //System.out.print(s);
                System.out.print(cplex.toString());
                cplex.output().println("Solution status = " + cplex.getStatus());
                System.out.println("h_d = " + cplex.getObjValue() + ", r=" + r);



				/*
				double[] val = cplex.getValues(x_dp[0]);
				//int ncols = cplex.getNcols();
				for (int i = 0; i <d_len; i++)
					for (int j = 0; j < p_len; ++j)
						cplex.output().println("D: " + i + " P:" + j + " Value = " + val[j]);

				for (int i = 0; i < d_len; i++)
					for (int j = 0; j < p_len; ++j) {
						System.out.println("V0: D: " + i + " P:" + j + " Value = " + delta_vdp[0][i][j]);
						System.out.println("V1: D: " + i + " P:" + j + " Value = " + delta_vdp[1][i][j]);
						//System.out.println("D: " + i + " P:" + j + " Value = " + delta_vdp[0][i][j]);
					}
					*/
            } else {
                System.out.println("Model not solved");
            }


        } catch (IloException exc) {
            exc.printStackTrace();
        }
    }


    private void fill_delta_s_t() {
        for (Integer v : V) {
            for (Integer d : D) {
                for (Path p : P) {
                    if (p.getD() == d) {
                        if (p.getStartV() == v) {
                            s_vd[v - 1][d - 1] = 0;
                        } else {
                            s_vd[v - 1][d - 1] = 1;
                        }
                        if (p.getEndV() == v) {
                            t_vd[v - 1][d - 1] = 0;
                        } else {
                            t_vd[v - 1][d - 1] = 1;
                        }
                        break;
                    }
                }
            }
        }
        //for (Integer v_p : p.getAllV()) {
        for (Integer d : D) {
            for (Integer v : V) {
                for (Path p : P) {
                    if (p.checkInPath(v)) {
                        if (p.getD() == d) {
                            delta_vdp[v - 1][d - 1][p.getNr() - 1] = 1;
                        } else {
                            delta_vdp[v - 1][d - 1][p.getNr() - 1] = 0;
                        }


                    } else {
                        delta_vdp[v - 1][d - 1][p.getNr() - 1] = 0;
                    }
                }
            }
        }
    }

    private void startData() {
        for (int i = 0; i <= v_len - 1; i++) {
            V[i] = i + 1;
        }
        for (int i = 0; i <= d_len - 1; i++) {
            D[i] = i + 1;
        }

        if (version == 3) {


            //Wersja V3
            P[0] = new Path(1, new Integer[]{1, 2});
            P[1] = new Path(1, new Integer[]{1, 3, 2});
            P[2] = new Path(2, new Integer[]{1, 3});
            P[3] = new Path(2, new Integer[]{1, 2, 3});
            P[4] = new Path(3, new Integer[]{2, 1});
            P[5] = new Path(3, new Integer[]{2, 3, 1});
            P[6] = new Path(4, new Integer[]{2, 3});
            P[7] = new Path(4, new Integer[]{2, 1, 3});
            P[8] = new Path(5, new Integer[]{3, 1});
            P[9] = new Path(5, new Integer[]{3, 2, 1});
            P[10] = new Path(6, new Integer[]{3, 2});
            P[11] = new Path(6, new Integer[]{3, 1, 2});

            T_v[0] = 80;
            T_v[1] = 100;
            T_v[2] = 120;

            R_v[0] = 30;
            R_v[1] = 30;
            R_v[2] = 60;

            pc_v[0] = 0.5;
            pc_v[1] = 0.6;
            pc_v[2] = 0.4;
        }
        if (version == 2) {
            //Wersja V2

            Parser p = new Parser();

            P[0] = new Path(1, p.compute(1, 2, 1));
            for (int i=0; i < P[0].V.length; i++) {
                System.out.println("P0: "+ P[0].V[i]);
            }

            P[1] = new Path(2, p.compute(2, 1, 1));
            for (int i=0; i < P[0].V.length; i++) {
                System.out.println("P0: "+ P[0].V[i]);
            }
            /*
            P[0] = new Path(1, new Integer[]{1, 2});
            P[1] = new Path(2, new Integer[]{2, 1});
*/
            T_v[0] = 30;
            T_v[1] = 40;


            R_v[0] = 50;
            R_v[1] = 60;


            pc_v[0] = 0.7;
            pc_v[1] = 0.8;
        }

        if (version == 4 ){
            Parser p = new Parser();
            /*
            for (int k=0; k<p_len; k++){
                for (int j=1; j<= v_len; j++){
                    for (int i=0; i<=v_len; i++){
                        P[k] = new Path(1, p.compute(, 2, 1));
                    }
                }
            }
            */
            P[0] = new Path(1, p.compute(1, 2, 1));
            P[1] = new Path(1, p.compute(1, 3, 1));
            P[2] = new Path(2, p.compute(1, 4, 1));
            P[3] = new Path(2, p.compute(2, 1, 1));
            P[4] = new Path(3, p.compute(2, 4, 1));
            P[5] = new Path(3, p.compute(3, 1, 1));
            P[6] = new Path(4, p.compute(3, 4, 1));
            P[7] = new Path(4, p.compute(4, 1, 1));
            P[8] = new Path(5, p.compute(4, 2, 1));
            P[9] = new Path(5, p.compute(4, 3, 1));
            P[10] = new Path(6, p.compute(2, 3, 1));
            P[11] = new Path(6, p.compute(3, 2, 1));
/*
            for (int j=0; j< p_len; j++) {
                for (int i = 0; i < P[1].getAllV().length; i++) {
                    System.out.println("P"+j+": " + P[j].V[i]);
                }
            }
*/
            T_v[0] = 70;
            T_v[1] = 60;
            T_v[2] = 60;
            T_v[3] = 40;

            R_v[0] = 60;
            R_v[1] = 90;
            R_v[2] = 40;
            R_v[3] = 50;

            pc_v[0] = 0.5;
            pc_v[1] = 0.8;
            pc_v[2] = 0.5;
            pc_v[3] = 0.4;
        }
    }

}
