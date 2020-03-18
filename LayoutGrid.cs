using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using MyExtension;

using UnityEngine.UI;

 public class LayoutGrid : MonoBehaviour
{

    [System.NonSerialized]
    public MeshRenderer[] MeshRenders;
    public Vector2[] Centers { get; }// 変更不可
    public Vector2 CentersMin { get; }
    public Vector2 CentersMax { get; }
    public int[] All_Index { get; }
    public int[] AliveIndex { get; set; }
    public int[][] AgentChild { get; set; }
    public eStatus[] Status { get; set; }
    public int[] Owner { get; set; } // -1 = this grid is not occupied eneyothers

    public float[] targetRatio { get; set; }
    public float[] currentRatio { get; set; }

    public float GridSize { get; } // m単位にする可能性があるためfloatが適切
    public int X_Extend { get; }
    public int Y_Extend { get; }
    public int AgentCount { get; }
    public int MaxIndex { get; }
    public int MinIndex { get; }
    public int CenterIndex { get; }
    public float Diagonal { get; }


    //IReadOnlyCollection
    //public IReadOnlyList<Entity> Items { get; }   // 変更不可
    public LayoutGrid(float gridSize, int gridNumX, int gridNumY, int agentCount)
    {

        this.GridSize = gridSize;
        this.X_Extend = gridNumX;
        this.Y_Extend = gridNumY;
        this.AgentCount = agentCount;
        this.AgentChild = new int[AgentCount][];
        for (int i = 0; i < AgentChild.Length; i++)
        {
            AgentChild[i] = new int[2];
        }
        this.Diagonal = gridSize * 1.5f; //root2 = 1.414

        this.MaxIndex = (gridNumX * gridNumY) -1;
        this.MinIndex = 0;
        this.CenterIndex = GetIndex(Mathf.FloorToInt(X_Extend / 2.0f), Mathf.FloorToInt(Y_Extend / 2.0f));

        MeshRenders = new MeshRenderer[MaxIndex+1];
        Centers = new Vector2[MaxIndex + 1];
        All_Index = new int[MaxIndex + 1];
        Status = new eStatus[MaxIndex + 1];
        Owner = new int[MaxIndex + 1];
        targetRatio = new float[AgentCount];
        currentRatio = new float[AgentCount];
        MakeTargetRatio(100f);

        int idx = 0;
        for (int i = 0; i < gridNumY; i++)
        {
            for (int j = 0; j < gridNumX; j++)
            {
                
                var center = new Vector2((gridSize * j) + (float)(gridSize / 2.0), (gridSize * i) + (float)(gridSize / 2.0));
                Centers[idx] = center;
                MeshRenders[idx] = MakeMesh(center, gridSize, idx.ToString()).GetComponent<MeshRenderer>();
                Status[idx] = eStatus.Alive;
                Owner[idx] = -1;
                All_Index[idx] = idx;
                idx++;
            }
        }
        CentersMin = Centers.First();
        CentersMax = Centers.Last();
    }

    /// <summary>
    /// Owner Statues をリセット
    /// </summary>
    public void ResetGrid()
    {
        for (int i = 0; i < Owner.Length; i++)
        {
            Owner[i] = -1;
            Status[i] = eStatus.Alive;
        }
    }

    /// <summary>
    /// scatterRate が1.1に近いと差がばらつきが大きなる。 大きな値だと均等になる
    /// </summary>
    /// <param name="scatterRate"></param>
    public void MakeTargetRatio(float scatterRate, float ratioSum = 1f)
    {
        if (scatterRate < 1f) scatterRate = 1f;
        if (ratioSum > 1f) ratioSum = 1f;

        float ratioAverage = ratioSum / (float)AgentCount;
        for (int i = 0; i < targetRatio.Length; i++)
        {
            targetRatio[i] = ratioAverage;
        }

        var rand = new System.Random();
        for (int i = 0; i < targetRatio.Length - 1; i++)
        {
            float randVal = (float)Remap(rand.NextDouble(), 0, 1.0, -ratioAverage / scatterRate, ratioAverage / scatterRate);
            targetRatio[i] += randVal;
            targetRatio[i + 1] -= randVal;
        }

        /*
        foreach (float val in targetRatio)
        {
            Debug.Log(val);
        }
        Debug.Log("sum" + targetRatio.Sum());
        */

    }
    /// <summary>
    ///  Agentの所有するグリッド数 / AliveIndex
    /// </summary>
    public void MakeCurrentRatio()
    {
        for (int i = 0; i < currentRatio.Length; i++)
        {
            int countOwnerI = Owner.Where(owner => owner == i).Count();

            currentRatio[i] = (float)countOwnerI / (float)AliveIndex.Length;
        }
    }


    public static double Remap(double value, double from1, double to1, double from2, double to2)
    {
        return (value - from1) / (to1 - from1) * (to2 - from2) + from2;
    }
    /// <summary>
    /// 中央からoffsetNum分離れたグリッドの中からランダムに選択し、コーナーとつなげてグリッドの形を変える。消されたらDeadになる。
    /// </summary>
    /// <param name="iteration"></param>
    /// <param name="offsetNum"></param>
    public void RandomReduce(int iteration, int offsetNum)
    {
        float range = offsetNum * Diagonal;
        List<int> maskedIdx = GetInRange(CenterIndex, range, All_Index.ToList());
        List<int> availableIdx = new List<int>();
        foreach (int idx in All_Index)
        {
            if (!maskedIdx.Contains(idx)) 
            {
                availableIdx.Add(idx);
            }
        }
        List<int> corners = GetCorners();

        for (int i = 0; i < iteration; i++)
        {
            availableIdx.Jitter();
            int randIdx = availableIdx[0];
            availableIdx.RemoveAt(0);

            var corner = ClosedGrid(randIdx, corners);
            var rectIdx = GetRect(randIdx, corner);

            foreach (int idx in rectIdx)
            {
                Status[idx] = eStatus.Dead;
            }
        }

        var aliveList = new List<int>();
        for (int j = 0; j < Status.Length; j++)
        {
            if (Status[j] == eStatus.Alive)
            {
                aliveList.Add(j);
            }
        }
        AliveIndex = aliveList.ToArray();
    }

    /// <summary>
    /// 二つのGridの距離を返す
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    private float GridDistance(int from, int to)
    {
        var vector = Centers[to] - Centers[from];
        return vector.magnitude;
    }

    /// <summary>
    /// form から　sCluster内の一番近いindexを返す
    /// </summary>
    /// <param name="from"></param>
    /// <param name="sCluster"></param>
    /// <returns></returns>
    public int ClosedGrid(int from, List<int> sCluster)
    {

        float[] distances = new float[sCluster.Count];
        int[] indexs = sCluster.ToArray();

        for (int i = 0; i < sCluster.Count; i++)
        {
            float dist = GridDistance(from, sCluster[i]);
            distances[i] = dist;
        }
        
        System.Array.Sort(distances, indexs);

        return indexs[0];

    }

    /// <summary>
    /// Grid全体の四隅のindexを得る
    /// </summary>
    /// <returns></returns>
    public List<int> GetCorners()
    {
        var rtnList = new List<int> {
            0,
            X_Extend - 1,
            MaxIndex,
           (Y_Extend - 1) * X_Extend
        };
        return rtnList;
    }

    /// <summary>
    /// range内にあるグリッドのindexを返す
    /// </summary>
    /// <param name="from"></param>
    /// <param name="range"></param>
    /// <returns></returns>
    private List<int> GetInRange(int from, float range, List<int> sCluster)
    {
        var rtnList = new List<int>() { from };

        for (int i = 0; i < sCluster.Count; i++)
        {
            var to = sCluster[i];
            if (GridDistance(from, to) < range)
            {
                rtnList.Add(to);
            }
        }

        return rtnList;
    }

    /// <summary>
    /// もしその方向が壁なら-1を返す
    /// </summary>
    /// <param name="current"></param>
    /// <returns></returns>
    public int GetRight(int current)
    {
       
        if (Centers[current].x == CentersMax.x) return -1; //一番右の列

        return current + 1;
    }

    /// <summary>
    /// もしその方向が壁なら-1を返す
    /// </summary>
    /// <param name="current"></param>
    /// <returns></returns>
    public int GetLeft(int current)
    {
        if (Centers[current].x == CentersMin.x) return -1; //一番左側

        return current - 1;
    }

    /// <summary>
    /// もしその方向が壁なら-1を返す
    /// </summary>
    /// <param name="current"></param>
    /// <returns></returns>
    public int GetUp(int current)
    {
        if (Centers[current].y == CentersMax.y) return -1; //一番左側

        return current + X_Extend;
    }

    /// <summary>
    /// もしその方向が壁なら-1を返す
    /// </summary>
    /// <param name="current"></param>
    /// <returns></returns>
    public int GetDown(int current)
    {
        if (Centers[current].y == CentersMin.y) return -1; //一番左側

        return current - X_Extend;
    }

    
    /// <summary>
    /// 右側になにがあるかと距離を、float少数（－1～1）二行1列で返す
    /// </summary>
    /// <param name="agentNum"></param>
    /// <returns></returns>
    public float[] GetRightOwner(int agentNum)
    {
        List<int> rightSide = new List<int>();
        float rightOwner = -1f;
        float rightOwnerDist = 0;

        //Childのどっちがｘの側の最大値をもっているのか判断してそこから右側の辺を得る
        Child xMaxChild = Centers[AgentChild[agentNum][(int)Child.A]].x > Centers[AgentChild[agentNum][(int)Child.B]].x ? Child.A : Child.B;
        Child xMinChild = xMaxChild == Child.A ? Child.B : Child.A;
        int maxYCount = GetY_Count(AgentChild[agentNum][(int)xMaxChild]);
        int minYCount = GetY_Count(AgentChild[agentNum][(int)xMinChild]);
       
        //右辺になるGridを追加していく
        int maxChildIdx = AgentChild[agentNum][(int)xMaxChild];
        rightSide.Add(maxChildIdx);
        //xMaxChildのほうがxMinChildよりｙ座標が大きい場合
        if (Centers[AgentChild[agentNum][(int)xMaxChild]].y > Centers[AgentChild[agentNum][(int)xMinChild]].y)
        {
            for (int i = 0; i < Mathf.Abs(maxYCount - minYCount); i++)
            {
                int down = GetDown(maxChildIdx);
                if (down == -1) break;
                rightSide.Add(down);
            }
        }
        else
        {
            for (int i = 0; i < Mathf.Abs(maxYCount - minYCount); i++)
            {
                int up = GetUp(maxChildIdx);
                if (up == -1) break;
                rightSide.Add(up);
            }
        }

        //rightSideないのGridを右に拡張していき、何回目に何があるか調べる
        for (int i = 0; i < X_Extend; i++)
        {
            for (int j = 0; j < rightSide.Count; j++)
            {
                int nextRight = GetRight(rightSide[j]);

                //もし領域外なら //もし壁にあたるなら
                if (nextRight == -1 || Status[nextRight] == eStatus.Dead)
                {
                    rightOwner = -1;
                    rightOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                //もしほかのエージェントの領域にあたるなら
                if (Owner[nextRight] != -1)
                {
                    //rightOwner = (float)Owner[nextRight] / (float)agentNum;
                    rightOwner = 1;
                    rightOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                rightSide[j] = nextRight;
            }
        }
    loopEnd: return new float[2] { rightOwner, rightOwnerDist };
    }

    public float[] GetLeftOwner(int agentNum)
    {
        List<int> leftSide = new List<int>();
        float leftOwner = -1f;
        float leftOwnerDist = 0;

        //Childのどっちがｘの側の最小値をもっているのか
        Child xMaxChild = Centers[AgentChild[agentNum][(int)Child.A]].x < Centers[AgentChild[agentNum][(int)Child.B]].x ? Child.A : Child.B;
        Child xMinChild = xMaxChild == Child.A ? Child.B : Child.A;
        int maxYCount = GetY_Count(AgentChild[agentNum][(int)xMaxChild]);
        int minYCount = GetY_Count(AgentChild[agentNum][(int)xMinChild]);

        //右辺になるGridを追加していく
        int maxChildIdx = AgentChild[agentNum][(int)xMaxChild];
        leftSide.Add(maxChildIdx);

        //xMaxChildのほうがxMinChildよりｙ座標が大きい場合
        if (Centers[AgentChild[agentNum][(int)xMaxChild]].y > Centers[AgentChild[agentNum][(int)xMinChild]].y)
        {
            for (int i = 0; i < Mathf.Abs(maxYCount - minYCount); i++)
            {
                int down = GetDown(maxChildIdx);
                if (down == -1) break;
                leftSide.Add(down);
            }
        }
        else
        {
            for (int i = 0; i < Mathf.Abs(maxYCount - minYCount); i++)
            {
                int up = GetUp(maxChildIdx);
                if (up == -1) break;
                leftSide.Add(up);
            }
        }


        //rightSideないのGridを右に拡張していき、何回目に何があるか調べる
        for (int i = 0; i < X_Extend; i++)
        {
            for (int j = 0; j < leftSide.Count; j++)
            {
                int nextRight = GetLeft(leftSide[j]);

                //もし領域外なら //もし壁にあたるなら
                if (nextRight == -1 || Status[nextRight] == eStatus.Dead)
                {
                    leftOwner = -1;
                    leftOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                //もしほかのエージェントの領域にあたるなら
                if (Owner[nextRight] != -1)
                {
                    //leftOwner = (float)Owner[nextRight] / (float)agentNum;
                    leftOwner = 1;
                    leftOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                leftSide[j] = nextRight;
            }
        }
    loopEnd: return new float[2] { leftOwner, leftOwnerDist };
    }

    public float[] GetUPOwner(int agentNum)
    {
        List<int> upSide = new List<int>();
        float upOwner = -1f;
        float upOwnerDist = 0;

        //Childのどっちがyの側の最大値をもっているのか
        Child yMaxChild = Centers[AgentChild[agentNum][(int)Child.A]].y > Centers[AgentChild[agentNum][(int)Child.B]].y ? Child.A : Child.B;
        Child yMinChild = yMaxChild == Child.A ? Child.B : Child.A;
        int maxXCount = GetX_Count(AgentChild[agentNum][(int)yMaxChild]);
        int minXCount = GetX_Count(AgentChild[agentNum][(int)yMinChild]);

        //右辺になるGridを追加していく
        int maxChildIdx = AgentChild[agentNum][(int)yMaxChild];
        upSide.Add(maxChildIdx);
        //
        if (Centers[AgentChild[agentNum][(int)yMaxChild]].x > Centers[AgentChild[agentNum][(int)yMinChild]].x)
        {
            for (int i = 0; i < Mathf.Abs(maxXCount - minXCount); i++)
            {
                int left = GetLeft(maxChildIdx);
                if (left == -1) break;
                upSide.Add(left);
            }
        }
        else
        {
            for (int i = 0; i < Mathf.Abs(maxXCount - minXCount); i++)
            {
                int right = GetRight(maxChildIdx);
                if (right == -1) break;
                upSide.Add(right);
            }
        }

        //rightSideないのGridを右に拡張していき、何回目に何があるか調べる
        for (int i = 0; i < X_Extend; i++)
        {
            for (int j = 0; j < upSide.Count; j++)
            {
                int next = GetUp(upSide[j]);

                //もし領域外なら //もし壁にあたるなら
                if (next == -1 || Status[next] == eStatus.Dead)
                {
                    upOwner = -1;
                    upOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                //もしほかのエージェントの領域にあたるなら
                if (Owner[next] != -1)
                {
                    //upOwner = (float)Owner[next] / (float)agentNum;
                    upOwner = 1;
                    upOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                upSide[j] = next;
            }
        }
    loopEnd: return new float[2] { upOwner, upOwnerDist };
    }

    public float[] GetDownOwner(int agentNum)
    {
        List<int> upSide = new List<int>();
        float upOwner = -1f;
        float upOwnerDist = 0;

        //Childのどっちがyの側の最大値をもっているのか
        Child yMaxChild = Centers[AgentChild[agentNum][(int)Child.A]].y < Centers[AgentChild[agentNum][(int)Child.B]].y ? Child.A : Child.B;
        Child yMinChild = yMaxChild == Child.A ? Child.B : Child.A;
        int maxXCount = GetX_Count(AgentChild[agentNum][(int)yMaxChild]);
        int minXCount = GetX_Count(AgentChild[agentNum][(int)yMinChild]);

        //右辺になるGridを追加していく
        int maxChildIdx = AgentChild[agentNum][(int)yMaxChild];
        upSide.Add(maxChildIdx);
        //
        if (Centers[AgentChild[agentNum][(int)yMaxChild]].x > Centers[AgentChild[agentNum][(int)yMinChild]].x)
        {
            for (int i = 0; i < Mathf.Abs(maxXCount - minXCount); i++)
            {
                int left = GetLeft(maxChildIdx);
                if (left == -1) break;
                upSide.Add(left);
            }
        }
        else
        {
            for (int i = 0; i < Mathf.Abs(maxXCount - minXCount); i++)
            {
                int right = GetRight(maxChildIdx);
                if (right == -1) break;
                upSide.Add(right);
            }
        }

        //rightSideないのGridを右に拡張していき、何回目に何があるか調べる
        for (int i = 0; i < X_Extend; i++)
        {
            for (int j = 0; j < upSide.Count; j++)
            {
                int next = GetDown(upSide[j]);

                //もし領域外なら //もし壁にあたるなら
                if (next == -1 || Status[next] == eStatus.Dead)
                {
                    upOwner = -1;
                    upOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                //もしほかのエージェントの領域にあたるなら
                if (Owner[next] != -1)
                {
                    //upOwner = (float)Owner[next] / (float)agentNum;
                    upOwner = 1;
                    upOwnerDist = (float)i / (float)X_Extend;
                    goto loopEnd;
                }

                upSide[j] = next;
            }
        }
    loopEnd: return new float[2] { upOwner, upOwnerDist };
    }

    /// <summary>
    /// 四方向のindexを取得する 
    /// </summary>
    /// <param name="current"></param>
    /// <returns></returns>
    public Dictionary<Axis, int> GetNeighbor(int current)
    {
        var neighbor = new Dictionary<Axis, int>
        {
            { Axis.Right, GetRight(current) },
            { Axis.Left, GetLeft(current) },
            { Axis.Down, GetDown(current) },
            { Axis.Up, GetUp(current) }
        };

        return neighbor;
    }
    /// <summary>
    ///生きてて壁じゃない方向を収集する。　所有者がだれかは無視。
    /// </summary>
    /// <param name="neighbor"></param>
    /// <returns></returns>
    public List<Axis> GetOpenNeighbor(int current)
    {
        List<Axis> openAxis = new List<Axis>();

        var neighborDict = GetNeighbor(current);
        int right = neighborDict[Axis.Right];
        int left = neighborDict[Axis.Left];
        int up = neighborDict[Axis.Up];
        int down = neighborDict[Axis.Down];

        if (right != -1 && Status[right] == eStatus.Alive)
        {
            openAxis.Add(Axis.Right);
        }

        if (left != -1 && Status[left] == eStatus.Alive)
        {
            openAxis.Add(Axis.Left);
        }

        if (up != -1 && Status[up] == eStatus.Alive)
        {
            openAxis.Add(Axis.Up);
        }

        if (down != -1 && Status[down] == eStatus.Alive)
        {
            openAxis.Add(Axis.Down);
        }

        return openAxis;
    }
    /// <summary>
    /// 2点間の四角形に入るgirdのインデックスを返す
    /// </summary>
    /// <param name="ptA"></param>
    /// <param name="ptB"></param>
    /// <returns></returns>
    public List<int> GetRect(int ptA, int ptB)
    {
        var rtnList = new List<int>();

       int agentA_X = GetX_Count(ptA);
       int agentA_Y = GetY_Count(ptA);
       int agentB_X = GetX_Count(ptB);
       int agentB_Y = GetY_Count(ptB);
        //Debug.Log(ptB + "ptb " + agentA_X + " " + agentA_Y + " " + agentB_X + " " + agentB_Y);

        int agentMaxX = agentA_X >= agentB_X ? agentA_X : agentB_X;
        int agentMaxY = agentA_Y >= agentB_Y ? agentA_Y : agentB_Y;

        int agentMinX = agentA_X <= agentB_X ? agentA_X : agentB_X;
        int agentMinY = agentA_Y <= agentB_Y ? agentA_Y : agentB_Y;

        //Debug.Log(ptB + "ptb " +agentMinX + " " + agentMinY + " " + agentMaxX + " " + agentMaxY);
        //intervalX = 0 ～ xGridNum
        int interValX = agentMaxX - agentMinX;
        int interValY = agentMaxY - agentMinY;

        for (int i = 0; i < interValY + 1; i++)
        {
            int y = agentMinY + i;
            for (int j = 0; j <  interValX + 1; j++)
            {
                int x = agentMinX + j;
         
                rtnList.Add(GetIndex(x, y));
            }
        }

        return rtnList;
    }

    public List<int> FirstPosition()
    {
        var initialPos = new List<int>();
        List<int> aliveBuff = AliveIndex.ToList();
     
        //適当な場所をエージェント分選ぶ
        for (int i = 0; i < 10000; i++)
        {
            aliveBuff.Jitter();
            //System.Array.Copy(aliveBuff, initialIdx, AgentCount);

            initialPos.Add(aliveBuff[0]);
            if (initialPos.Count == AgentCount)
            {
                break;
            }
            aliveBuff.RemoveAt(0);


            // List<int> inRange = GetInRange(initialPos.Last(), minDist, aliveBuff) ;
            var neighbor = GetNeighbor(initialPos.Last());

            foreach (var idx in neighbor.Values)
            {
                if (aliveBuff.Contains(idx))
                {
                    aliveBuff.Remove(idx);
                }
            }
            if (aliveBuff.Count == 0)
            {
                initialPos.Clear();
                aliveBuff = AliveIndex.ToList();
                Debug.Log("FirstPosition iteration Reset");
            }
        }

        return initialPos;
    }


    public void SplitPosition(List<int> firstPosition)
    { 
        for(int i = 0; i < firstPosition.Count; i++)
        {
            var neighbor = GetNeighbor(firstPosition[i]);
            int right = neighbor[Axis.Right];
            int left = neighbor[Axis.Left];
            int up = neighbor[Axis.Up];
            int down = neighbor[Axis.Down];

            Owner[firstPosition[i]] = i;
            AgentChild[i][(int)Child.A] = firstPosition[i];

            if (right != -1 && Status[right] != eStatus.Dead && Owner[right] == -1)
            {
                Owner[right] = i;
                AgentChild[i][(int)Child.B] = right;

                continue;
            }

            if (left != -1 && Status[left] != eStatus.Dead && Owner[left] == -1)
            {
                Owner[left] = i;
                AgentChild[i][(int)Child.B] = left;
                continue;
            }

            if (up != -1 && Status[up] != eStatus.Dead && Owner[up] == -1)
            {
                Owner[up] = i;
                AgentChild[i][(int)Child.B] = up;
                continue;
            }

            if (down != -1 && Status[down] != eStatus.Dead && Owner[down] == -1)
            {
                Owner[down] = i;
                AgentChild[i][(int)Child.B] = down;
                continue;
            }
            AgentChild[i][(int)Child.B] = AgentChild[i][(int)Child.A];
        }
        return;
    }

    /// <summary>
    /// 行きたい方向に移動させてステータスを変える。
    /// もしほかのエージェントの領域ならオフセット可能なら自動でしてくれる。できなければ移動しない。
    /// </summary>
    /// <param name="agentNum"></param>
    /// <param name="current"></param>
    /// <param name="next"></param>
    public void AgentMove(int agentNum, Child child, Axis nextAxis)
    {

        Child anotherChild = child == Child.A ? Child.B : Child.A;
        if (child == Child.A)
        {
            anotherChild = Child.B;
        }
        else
        {
            anotherChild = Child.A;
        }

        int anotherChildIdx = AgentChild[agentNum][(int)anotherChild];

        int current = AgentChild[agentNum][(int)child];
        int next = -1;

        switch (nextAxis)
        {
            case Axis.None:
                return;

            case Axis.Right:
                next = GetRight(current);
                break;

            case Axis.Left:
                next = GetLeft(current);
                break;

            case Axis.Up:
                next = GetUp(current);
                break;

            case Axis.Down:
                next = GetDown(current);
                break;
        }
        if (next == -1 || Status[next] == eStatus.Dead ) return; // 進行方向が壁 か死んでる

        //現在のエージェントの拡張されるマスの中にほかのエージェントがいないかどうかを調べる。ついでに次のマスの所有者も知っておく
        int preAgentNum = -1;
        int count = 0;
        List<int> nextRect = GetRect(anotherChildIdx, next);
        for (int i = 0; i < nextRect.Count; i++)
        {
            var idx = nextRect[i];
            //生きててほかのエージェントの支配下なら
            if (Owner[idx] != agentNum && Owner[idx] != -1 && Status[idx] == eStatus.Alive)
            {
                preAgentNum = Owner[idx];
                count++;
            }
        }
        if (count > 1) return;

        //空いてるなら普通に移動する。
        if (preAgentNum == -1 && count == 0)
        {
            //Childの値を更新
            AgentChild[agentNum][(int)child] = next;

            //現在のGridのOwnerを-1にして次のオーナーになる。
            for (int i = 0; i < Owner.Length; i++)
            {
                if (Owner[i] == agentNum)
                {
                    Owner[i] = -1;
                }
            }

            //範囲ないにあるGridのオーナーになる
            List<int> inRect = GetRect(AgentChild[agentNum][0], AgentChild[agentNum][1]);
            foreach (int idx in inRect)
            {
                if (Status[idx] != eStatus.Dead)
                {
                    Owner[idx] = agentNum;
                }
            }
            return;
        }


        //進む方向のオーナーがラインならなにもしない。
        if (isLine(preAgentNum))
        {
            return;
        }
        //次すすmところにいるownerがlineではない
        else
        {
            //進む先のOwnerをオフセットする。
            int preOwnerIdxA = AgentChild[preAgentNum][(int)Child.A];
            int preOwnerIdxB = AgentChild[preAgentNum][(int)Child.B];
            int preOwnerNextA = -1;
            int preOwnerNextB = -1;

            switch (nextAxis)
            {
                case Axis.None:
                    return;

                case Axis.Right:
                    preOwnerNextA = GetRight(preOwnerIdxA);
                    preOwnerNextB = GetRight(preOwnerIdxB);
                    break;

                case Axis.Left:
                    preOwnerNextA = GetLeft(preOwnerIdxA);
                    preOwnerNextB = GetLeft(preOwnerIdxB);
                    break;

                case Axis.Up:
                    preOwnerNextA = GetUp(preOwnerIdxA);
                    preOwnerNextB = GetUp(preOwnerIdxB);
                    break;

                case Axis.Down:
                    preOwnerNextA = GetDown(preOwnerIdxA);
                    preOwnerNextB = GetDown(preOwnerIdxB);
                    break;
            }

            //領域外でも死んでもないなら
            if (preOwnerNextA != -1 && preOwnerNextB != -1 && 
                Status[preOwnerNextA] == eStatus.Alive && Status[preOwnerNextB] == eStatus.Alive)
            {
                var preOwnerNewRect = GetRect(preOwnerNextA, preOwnerNextB);
                foreach (int idx in preOwnerNewRect)
                {
                    //次のグリッドが誰かのものなら
                    if (Owner[idx] != preAgentNum && Owner[idx] != -1)
                    {
                        return;
                    }
                }
            }
            else
            { 
                return; 
            }

            //現在のAgentの次の領域内にプレえじぇんと以外がいたらだめ,だけどこれは最初にやってる

            //前任者のAgentを動かす
            AgentChild[preAgentNum][(int)Child.A] = preOwnerNextA;
            AgentChild[preAgentNum][(int)Child.B] = preOwnerNextB;

            //現在のAgentを動かす。
            AgentChild[agentNum][(int)child] = next;

            //現在のGridのOwnerを-1にして次のオーナーになる。
            for (int i = 0; i < Owner.Length; i++)
            {
                if (Owner[i] == agentNum || Owner[i] == preAgentNum)
                {
                    Owner[i] = -1;
                }
            }
         
            //Rectを二匹更新する。
            List<int> inRect = GetRect(AgentChild[agentNum][0], AgentChild[agentNum][1]);
            foreach (int idx in inRect)
            {
                if (Status[idx] == eStatus.Alive)
                {
                    Owner[idx] = agentNum;
                }
            }

            inRect = GetRect(AgentChild[preAgentNum][0], AgentChild[preAgentNum][1]);
            foreach (int idx in inRect)
            {
                if (Status[idx] == eStatus.Alive)
                {
                    Owner[idx] = preAgentNum;
                }
            }
        }

        return;
    }

    private bool isLine(int agentNum)
    {
        if (agentNum == -1)
        {
            Debug.Log("isLine rececived -1");
            return false;
        }
        int childA = AgentChild[agentNum][0];
        int childB = AgentChild[agentNum][1];

        //agent child が同じ行か列にある。
        if (Centers[childA].x == Centers[childB].x || Centers[childA].y == Centers[childB].y)
        {
            return true;
        }

        return false;
    }

    private void SetAgent(int agentNumber, int gridIndex, int child)
    {
        //Status[gridIndex] = eStatus.Closed;
        Owner[gridIndex] = agentNumber;
        AgentChild[agentNumber][child] = gridIndex;
        return;
    }

    private Axis GetAxis(int from, int to)
    {
        Vector2 direction = (Centers[to] - Centers[from]).normalized;
        
        float dot = Vector2.Dot(direction, Vector2.right);
     
        //right
        if (1 > dot && dot > 0.7)
        {
            return Axis.Right;
        }

        //left
        if (-0.7 > dot && dot > -1)
        {
            return Axis.Left;
        }

        if (direction.y >= 0)
        {
            return Axis.Up;
        }

        return Axis.Down;
    }

    /// <summary>
    /// /grid の idx番目が　x軸に何個目か（左から数えて何行目か　（1～　xExtend)
    /// </summary>
    /// <param name="idx"></param>
    /// <returns></returns>
    private int GetX_Count(int idx)
    {
        if (idx == 0) return 1;

        int quotient = Mathf.FloorToInt(idx / X_Extend);
        return idx - (quotient * X_Extend) + 1;
    }

    /// <summary>
    /// /grid の idx番目のｙ座標（ｙ下から数えて何列目か　（1～　xExtend)
    /// </summary>
    /// <param name="idx"></param>
    /// <returns></returns>
    private int GetY_Count(int idx)
    {
        if (idx == 0) return 1;

        int quotient = Mathf.FloorToInt(idx / X_Extend);
        return quotient + 1;
    }

   
    /// <summary>
    /// (xCount, yCount)　の位置にあるグリッドのindexを返す
    /// </summary>
    /// <param name="xCount"></param>
    /// <param name="yCount"></param>
    /// <returns></returns>
    private int GetIndex(int xCount, int yCount)
    {
        return X_Extend * (yCount - 1) + xCount -1;
    }

 
    private GameObject MakeMesh(Vector2 center, float gridSize, string name = "")
    {
        GameObject meshObj = new GameObject("Mesh" + name);
        var mesh = new Mesh();
        
        mesh.vertices = new Vector3[] {
            new Vector3( center.x - (float)(gridSize /2.0), center.y - (float)(gridSize /2.0), 0),
            new Vector3(center.x - (float)(gridSize / 2.0), center.y + (float)(gridSize / 2.0), 0),
            new Vector3( center.x + (float)(gridSize /2.0), center.y - (float)(gridSize /2.0), 0),
            new Vector3( center.x + (float)(gridSize /2.0), center.y + (float)(gridSize /2.0), 0),
        };

        mesh.triangles = new int[] {
        0, 1, 2,
        1, 3, 2,
         };

        var meshRender = meshObj.AddComponent<MeshRenderer>();
        var meshFilter = meshObj.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;
        meshRender.material.color = Color.clear;
        meshRender.material.shader = Shader.Find("UI/Default");
    
        return meshObj;
    }


    public void PreviewStatues()
    {
       
        for (int i = 0; i < Status.Length; i++)
        {
            if (Status[i] == eStatus.Dead)
            {
                MeshRenders[i].material.color = Color.black;
            }
            else if(Status[i] == eStatus.Alive)
            {
                MeshRenders[i].material.color = Color.white;
            }
        }
    }

    public void PreviewOwner(Color[] colors)
    {
        for (int i = 0; i < Owner.Length; i++)
        {
            if (Owner[i] == -1 && Status[i] == eStatus.Alive)
            {
                MeshRenders[i].material.color = Color.white;
            }
            else if (Owner[i] != -1 && Status[i] == eStatus.Alive)
            {
                MeshRenders[i].material.color = colors[Owner[i]];
            }
          
            /*MeshRenders[i].transform.position = new Vector3(0, 0, 0);
            MeshRenders[i].transform.name = Owner[i].ToString() +" " + i.ToString() ;*/
        }
      

        /*
        for (int i = 0; i < AgentChild.Length; i++)
        {
            MeshRenders[AgentChild[i][(int)Child.A]].transform.position = new Vector3(0, 0, 1f);
            MeshRenders[AgentChild[i][(int)Child.B]].transform.position = new Vector3(0, 0, 1f);
        }
         */
    }
}

