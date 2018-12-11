import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

import java.util.ArrayList;
import java.util.Scanner;



public class Solver {
	
	
	
	public static double assignmentWithOrWithoutConstraintsHDB(boolean capacityConstraints){
		
			
		ArrayList<GRBVar[][]> vars_x = new ArrayList<GRBVar[][]>();
		for (int b=0;b<MainTests.numBlocks;b++){
			vars_x.add(new GRBVar[MainTests.numAgents][MainTests.numAptPerBlock[b]]);
		}
		try {

			//Model Creation :
			GRBEnv env = new GRBEnv();
			GRBModel model = new GRBModel(env);

			//Variables Creation :
			if (capacityConstraints){
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int i=0;i<MainTests.numAgents;i++){
						double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
						for (int j=0;j<numApt_b;j++){
							vars_x_b[i][j] = model.addVar(0, 1, utilities_b[j], GRB.BINARY, "xij");
						}
					}
				}
			}else{//Without capacity constraints, the optimal solution can be obtained with continuous variables:
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int i=0;i<MainTests.numAgents;i++){
						double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
						for (int j=0;j<numApt_b;j++){
							vars_x_b[i][j] = model.addVar(0, 1, utilities_b[j], GRB.CONTINUOUS, "xij");
						}
					}
				}
			}
			model.update();

			if (capacityConstraints){
				//Type-block constraints:
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int[] capacities_b = MainTests.capacities.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					int cpt = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						int numAgents_l = MainTests.numAgentPerType[l];
						GRBLinExpr expr = new GRBLinExpr();
						for (int i=0;i<numAgents_l;i++){
							GRBVar[] vars_x_b_i = vars_x_b[cpt];
							for (int j=0;j<numApt_b;j++){
								expr.addTerm(1, vars_x_b_i[j]);
							}
							cpt = cpt + 1;
						}
						model.addConstr(expr,GRB.LESS_EQUAL,capacities_b[l],"capacities");
					}
					model.update();	
				}
			}

			//matching constraints:
			for (int i=0;i<MainTests.numAgents;i++){
				GRBLinExpr expr = new GRBLinExpr();
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[] vars_x_b_i = vars_x.get(b)[i];
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int j=0;j<numApt_b;j++){
						expr.addTerm(1, vars_x_b_i[j]);
					}
				}
				model.addConstr(expr,GRB.LESS_EQUAL,1,"sortie");
			}
			model.update();
			for (int b=0;b<MainTests.numBlocks;b++){
				int numApt_b = MainTests.numAptPerBlock[b];
				GRBVar[][] vars_x_b = vars_x.get(b);
				for (int j=0;j<numApt_b;j++){
					GRBLinExpr expr = new GRBLinExpr();
					for (int i=0;i<MainTests.numAgents;i++){
						expr.addTerm(1, vars_x_b[i][j]);
					}
					model.addConstr(expr,GRB.LESS_EQUAL,1,"entree");
				}	
			}
			model.update();
		
			//Optimize:
			model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
			model.getEnv().set(GRB.IntParam.OutputFlag, 0);
			model.update();
			model.optimize();
			
			int optimstatus = model.get(GRB.IntAttr.Status);

		    if (optimstatus == GRB.Status.INFEASIBLE){
		    	model.dispose();
				env.dispose();
		    	return -Double.MAX_VALUE;
		    }else{
		    	double valSolution = model.get(GRB.DoubleAttr.ObjVal);  			    
		    	model.dispose();
				env.dispose();
				return valSolution;
		    }
		} catch (GRBException e) {
			//e.printStackTrace();
			return -Double.MAX_VALUE;
		}
	}
	
	
	
	
	public static double assignmentWithOrWithoutConstraintsChicago(boolean capacityConstraints){
		
		ArrayList<GRBVar[][]> vars_x = new ArrayList<GRBVar[][]>();
		for (int b=0;b<MainTests.numBlocks;b++){
			vars_x.add(new GRBVar[MainTests.numAgents][MainTests.numAptPerBlock[b]]);
		}
		try {
			//Model Creation:
			GRBEnv env = new GRBEnv();
			GRBModel model = new GRBModel(env);
			
			//Variables Creation:
			if (capacityConstraints){
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int i=0;i<MainTests.numAgents;i++){
						double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
						for (int j=0;j<numApt_b;j++){
							vars_x_b[i][j] = model.addVar(0, 1, utilities_b[j], GRB.BINARY, "xij");
						}
					}
				}
			}else{//Without capacity constraints, the optimal solution can be obtained with continuous variables:
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int i=0;i<MainTests.numAgents;i++){
						double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
						for (int j=0;j<numApt_b;j++){
							vars_x_b[i][j] = model.addVar(0, 1, utilities_b[j], GRB.CONTINUOUS, "xij");
						}
					}
				}
			}
			model.update();
			
			if (capacityConstraints){
				//Type-block constraints:
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[][] vars_x_b = vars_x.get(b);
					int[] capacities_b = MainTests.capacities.get(b);
					int numApt_b = MainTests.numAptPerBlock[b];
					int cpt = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						int numAgents_l = MainTests.numAgentPerType[l];
						GRBLinExpr expr = new GRBLinExpr();
						for (int i=0;i<numAgents_l;i++){
							GRBVar[] vars_x_b_i = vars_x_b[cpt];
							for (int j=0;j<numApt_b;j++){
								expr.addTerm(1, vars_x_b_i[j]);
							}
							cpt = cpt + 1;
						}
						model.addConstr(expr,GRB.LESS_EQUAL,capacities_b[l],"capacities");
					}
					model.update();	
				}
			}

			//matching constraints:
			for (int i=0;i<MainTests.numAgents;i++){
				GRBLinExpr expr = new GRBLinExpr();
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[] vars_x_b_i = vars_x.get(b)[i];
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int j=0;j<numApt_b;j++){
						expr.addTerm(1, vars_x_b_i[j]);
					}
				}
				model.addConstr(expr,GRB.LESS_EQUAL,1,"sortie");
			}
			model.update();
			for (int b=0;b<MainTests.numBlocks;b++){
				int numApt_b = MainTests.numAptPerBlock[b];
				GRBVar[][] vars_x_b = vars_x.get(b);
				for (int j=0;j<numApt_b;j++){
					GRBLinExpr expr = new GRBLinExpr();
					for (int i=0;i<MainTests.numAgents;i++){
						expr.addTerm(1, vars_x_b[i][j]);
					}
					model.addConstr(expr,GRB.LESS_EQUAL,1,"entree");
				}	
			}
			model.update();
			
			//Optimize:
			model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
			model.getEnv().set(GRB.IntParam.OutputFlag, 0);
			model.update();
			model.optimize();
			
			int optimstatus = model.get(GRB.IntAttr.Status);

		    if (optimstatus == GRB.Status.INFEASIBLE){	
		    	model.dispose();
				env.dispose();
				System.out.println("No solution in affectation.");
				Scanner sc = new Scanner(System.in);
				String str = sc.nextLine();
		    	return -Double.MAX_VALUE;
		    }else{   	
		    	double valSolution = model.get(GRB.DoubleAttr.ObjVal);   			    
		    	model.dispose();
				env.dispose();
				return valSolution;
		    }
		} catch (GRBException e) {
			//e.printStackTrace();
			return -Double.MAX_VALUE;
		}
	}
	
	
	
	
	
	public static double computeDeltaStar(double optVal){
		
		ArrayList<GRBVar[][]> vars_x = new ArrayList<GRBVar[][]>();
		for (int b=0;b<MainTests.numBlocks;b++){
			vars_x.add(new GRBVar[MainTests.numAgents][MainTests.numAptPerBlock[b]]);
		}
		try {
			//Model Creation:
			GRBEnv env = new GRBEnv();
			GRBModel model = new GRBModel(env);
			
			//Variables Creation:
			GRBVar t = model.addVar(-Double.MAX_VALUE, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "t");			
			for (int b=0;b<MainTests.numBlocks;b++){
				GRBVar[][] vars_x_b = vars_x.get(b);
				int numApt_b = MainTests.numAptPerBlock[b];
				for (int i=0;i<MainTests.numAgents;i++){
					for (int j=0;j<numApt_b;j++){
						vars_x_b[i][j] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "xij");
					}
				}
			}
			model.update();
			
			//We want the matching to have a utility equal to valOpt:
			GRBLinExpr expr = new GRBLinExpr();
			for (int b=0;b<MainTests.numBlocks;b++){
				GRBVar[][] vars_x_b = vars_x.get(b);
				int numApt_b = MainTests.numAptPerBlock[b];
				for (int i=0;i<MainTests.numAgents;i++){
					double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
					for (int j=0;j<numApt_b;j++){
						if (utilities_b[j]!=0)
							expr.addTerm(utilities_b[j]*Math.pow(10,5), vars_x_b[i][j]);
					}
				}
			}
			model.addConstr(expr,GRB.GREATER_EQUAL,(optVal-Math.pow(10, -10))*Math.pow(10,5),"opt");
			model.update();	
			
			//Min linearization:
			int cpt =0;
			GRBLinExpr exprt = new GRBLinExpr();
			exprt.addTerm(optVal*Math.pow(10,8)*1.0/MainTests.numAgents, t);
			for (int l=0;l<MainTests.numTypes;l++){
				expr = new GRBLinExpr();
				int numAgents_l = MainTests.numAgentPerType[l];
				int start = cpt;
				for (int i=start;i<start+numAgents_l;i++){
					for (int b=0;b<MainTests.numBlocks;b++){
						GRBVar[][] vars_x_b = vars_x.get(b);
						double[] utilities_b = MainTests.agents.get(i).getUtilities().get(b);
						int numApt_b = MainTests.numAptPerBlock[b];
						for (int j=0;j<numApt_b;j++){
							if (utilities_b[j]!=0)
								expr.addTerm(utilities_b[j]*Math.pow(10,8)*1.0/numAgents_l, vars_x_b[i][j]);
						}
					}
					cpt = cpt+1;
				}
				model.addConstr(exprt,GRB.LESS_EQUAL,expr,"min");
			}
			model.update();

			//matching constraints:
			for (int i=0;i<MainTests.numAgents;i++){
				expr = new GRBLinExpr();
				for (int b=0;b<MainTests.numBlocks;b++){
					GRBVar[] vars_x_b_i = vars_x.get(b)[i];
					int numApt_b = MainTests.numAptPerBlock[b];
					for (int j=0;j<numApt_b;j++){
						expr.addTerm(1, vars_x_b_i[j]);
					}
				}
				model.addConstr(expr,GRB.LESS_EQUAL,1,"sortie");
			}
			model.update();
			for (int b=0;b<MainTests.numBlocks;b++){
				int numApt_b = MainTests.numAptPerBlock[b];
				GRBVar[][] vars_x_b = vars_x.get(b);
				for (int j=0;j<numApt_b;j++){
					expr = new GRBLinExpr();
					for (int i=0;i<MainTests.numAgents;i++){
						expr.addTerm(1, vars_x_b[i][j]);
					}
					model.addConstr(expr,GRB.LESS_EQUAL,1,"entree");
				}	
			}
			model.update();
			
					
			//Optimize:
			model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
			model.getEnv().set(GRB.IntParam.OutputFlag, 0);
			model.update();
			model.optimize();
			
			int optimstatus = model.get(GRB.IntAttr.Status);

		    if (optimstatus == GRB.Status.INFEASIBLE){
		    	model.dispose();
				env.dispose();
		    	return -Double.MAX_VALUE;
		    }else{
		    	
		    	double valSolution = model.get(GRB.DoubleAttr.ObjVal);	
		    	/*for (int b=0;b<MainTests.numBlocks;b++){
					int numApt_b = MainTests.numAptPerBlock[b];
					double[][] solution = model.get(GRB.DoubleAttr.X,vars_x.get(b));
					for (int i=0;i<MainTests.numAgents;i++){
						for (int j=0;j<numApt_b;j++){
							if (solution[i][j]!= 0 && solution[i][j]!= 1){
					    		System.out.println(solution[i][j]);
							}
						}
					}
		    	}*/
		    	model.dispose();
				env.dispose();
				return valSolution;
		    }
		} catch (GRBException e) {
			//e.printStackTrace();
			return -Double.MAX_VALUE;
		}
	}
	

}