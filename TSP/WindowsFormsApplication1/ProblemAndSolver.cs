using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        private class Node : IComparable<Node>
        {
            private double bssf;
            private double[,] matrix;
            private ArrayList partialRoute;
            private int num;

            public Node()
            {
                bssf = 0;
                matrix = new double[0, 0];
                partialRoute = new ArrayList();
                num = -1;
            }

            public Node(double bssf, double[,] matrix, ArrayList partialRoute)
            {
                this.bssf = bssf;
                this.matrix = matrix;
                this.partialRoute = new ArrayList(partialRoute);
            }

            public double getbssf()
            {
                return bssf;
            }

            public double[,] getMatrix()
            {
                return matrix;
            }

            public ArrayList getRoute()
            {
                return partialRoute;
            }

            public void setbssf(double bssf)
            {
                this.bssf = bssf;
            }

            public void setMatrix(double[,] matrix)
            {
                this.matrix = matrix;
            }

            public void setRoute(ArrayList partialRoute)
            {
                this.partialRoute = partialRoute;
            }

            public int getNum()
            {
                return num;
            }

            public void setNum(int num)
            {
                this.num = num;
            }

            public int CompareTo(Node two)
            {
                //if ((this.getbssf() / (double)Cities.Length) < (two.getbssf() / (double)Cities.Length))
                //  return -1;
                //return 1;
                if (this.getRoute().Count > two.getRoute().Count)
                    return -1;
                else if (this.getRoute().Count == two.getRoute().Count)
                {
                    if (this.getbssf() <= two.getbssf())
                        return -1;
                    else
                        return 1;
                }
                else
                    return 1;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        private Node reduceMatrix(Node node)
        {
            double[,] updated = new double[Cities.Length, Cities.Length];
            double updatedbssf = node.getbssf();
            for(int i = 0; i < Cities.Length; i++)
            {
                for (int k = 0; k < Cities.Length; k++)
                    updated[i, k] = node.getMatrix()[i, k];
            }

            for(int i = 0; i < Cities.Length; i++)     //for each row
            {
                double min = double.PositiveInfinity;
                for(int k = 0; k < Cities.Length; k++)     //get minimum for each row
                {
                    if (updated[i, k] < min)
                        min = updated[i, k];
                }
                if(min != 0 && min != double.PositiveInfinity)
                {
                    for(int k = 0; k < Cities.Length; k++) //subtract minimum from each value
                        updated[i, k] -= min;
                    updatedbssf += min;
                }
            }

            for (int i = 0; i < Cities.Length; i++)     //for each column
            {
                double min = double.PositiveInfinity;
                for (int k = 0; k < Cities.Length; k++)     //get minimum for each column
                {
                    if (updated[k, i] < min)
                        min = updated[k, i];
                }
                if (min != 0 && min != double.PositiveInfinity)
                {
                    for (int k = 0; k < Cities.Length; k++) //subtract minimum from each value
                        updated[k, i] -= min;
                    updatedbssf += min;
                }
            }
            Node newNode = new Node(updatedbssf, updated, node.getRoute());
            newNode.setNum(node.getNum());
            return newNode;
        }

        private Node getPartialRoute(Node node, int rowNum, int colNum)
        {
            double[,] updated = new double[Cities.Length, Cities.Length];
            for(int i = 0; i < Cities.Length; i++)
            {
                for (int k = 0; k < Cities.Length; k++)
                    updated[i, k] = node.getMatrix()[i, k];
            }
            ArrayList partialRoute = new ArrayList();
            for (int i = 0; i < node.getRoute().Count; i++)
                partialRoute.Add(node.getRoute()[i]);

            double newbssf = node.getbssf() + node.getMatrix()[rowNum, colNum];
            for(int i = 0; i < Cities.Length; i++)
            {
                updated[rowNum, i] = double.PositiveInfinity;
                updated[i, colNum] = double.PositiveInfinity;
            }
            updated[rowNum, colNum] = double.PositiveInfinity;
            updated[colNum, rowNum] = double.PositiveInfinity;
            partialRoute.Add(Cities[colNum]);
            Node newNode = new Node(newbssf, updated, partialRoute);
            newNode.setNum(colNum);
            return reduceMatrix(newNode);
        }
        
        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];
            string[] origResults = defaultSolveProblem();
            double origbssf = Convert.ToDouble(origResults[COST]);
            double[,] originalDistances = new double[Cities.Length, Cities.Length];
            List<Node> todoNodes = new List<Node>();
            int numSolutions = 0;
            int numChildStates = 0;
            int statesPruned = 0;
            int numBSSF = 0;
            int maxNumOfStates = 0;
            double bestSolution = origbssf;
            ArrayList routeSolution = new ArrayList();

            for (int i = 0; i < Cities.Length; i++) //make original matrix from default route
            {
                for(int k = 0; k < Cities.Length; k++)
                {
                    if (i == k)
                        originalDistances[i, k] = double.PositiveInfinity;
                    else
                        originalDistances[i, k] = Cities[i].costToGetTo(Cities[k]);
                }
            }

            Stopwatch timer = new Stopwatch();
            timer.Start();
            routeSolution.Add(Cities[0]);
            Node startNode = new Node(0, originalDistances, routeSolution);
            //put nodes in queue
            todoNodes.Add(reduceMatrix(startNode));

            while(todoNodes.Count > 0)  //while the list is not empty
            {
                if (timer.Elapsed.TotalSeconds >= TIME_LIMIT)
                    break;
                Node node = todoNodes[0]; //get first element
                if(node.getRoute().Count == Cities.Length)  //if a leaf node
                {
                    numSolutions++;
                    TSPSolution sol = new TSPSolution(node.getRoute());
                    if (sol.costOfRoute() < bestSolution)   //reset the best solution if it is lower
                    {
                        bestSolution = sol.costOfRoute();
                        routeSolution = node.getRoute();
                        numBSSF++;
                    }
                }
                else    //if not a leaf node 
                {
                    if (node.getbssf() <= bestSolution)  //pruning not needed
                    {
                        for (int i = 1; i < Cities.Length; i++)
                        {
                            if (node.getMatrix()[node.getNum(), i] != double.PositiveInfinity)
                            {
                                todoNodes.Add(getPartialRoute(node, node.getNum(), i));
                                numChildStates++;
                            }
                        }
                    }
                    else
                        statesPruned++;
                }
                
                if (todoNodes.Count > maxNumOfStates)
                    maxNumOfStates = todoNodes.Count;
                todoNodes.RemoveAt(0);  //remove the node when finished calculating
                todoNodes.Sort();   //sorting based on priority
            }
            timer.Stop();
            bssf = new TSPSolution(routeSolution);
            Console.WriteLine("Number of child states is " + numChildStates);
            Console.WriteLine("Number of states pruned is " + statesPruned);
            Console.WriteLine("Number of bssf updates is " + numBSSF);
            Console.WriteLine("Number of stored states is " + maxNumOfStates);

            results[COST] = bssf.costOfRoute().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = numSolutions.ToString();
            return results;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for a greedy solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];

            // implement simulated annealing

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }
        #endregion
    }

}
