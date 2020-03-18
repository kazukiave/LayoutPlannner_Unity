using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MLAgents;
using MyExtension;
using System.Linq;
using MLAgents.Sensor;
using System.Text;
public class LayoutAgent : Agent
{
    [Header("set vectorGrid on same hieralchy")]
    private LayoutEnviroment layoutEnviroment;

    [System.NonSerialized]
    private LayoutGrid layoutGrid;

    private int agentNum;


    public bool showDebug = false;
    private StringBuilder obsVectors = new StringBuilder();

    public override void InitializeAgent()
    {
     
        layoutEnviroment = transform.parent.GetComponent<LayoutEnviroment>();
        layoutEnviroment.InitEnviroment();

        layoutGrid = layoutEnviroment.layoutGrid;
        agentNum = int.Parse(gameObject.name);
    }

    public override void CollectObservations()
    {
        int childA = layoutGrid.AgentChild[agentNum][(int)Child.A];
        int childB = layoutGrid.AgentChild[agentNum][(int)Child.B];
  
        Vector2 maxVector = new Vector2((float)layoutGrid.X_Extend, (float)layoutGrid.Y_Extend);
        Vector2 direction = (layoutGrid.Centers[childA] - layoutGrid.Centers[childB]) / maxVector;

        AddVectorObs(direction);
     
        //他のエージェントが所有している場所との距離と方向
        float[] rightObs = layoutGrid.GetRightOwner(agentNum);
        float[] leftObs = layoutGrid.GetLeftOwner(agentNum);
        float[] upObs = layoutGrid.GetUPOwner(agentNum);
        float[] downObs = layoutGrid.GetDownOwner(agentNum);
        AddVectorObs(rightObs);
        AddVectorObs(leftObs);
        AddVectorObs(upObs);
        AddVectorObs(downObs);

        //壁までの距離と方向

        //現在の面積と目的の面積の差分
        layoutGrid.MakeCurrentRatio();
        float ratioDiff = layoutGrid.targetRatio[agentNum] - layoutGrid.currentRatio[agentNum];//この差分が小さいほど評価が高い
    
        AddVectorObs(ratioDiff);
    

        if (showDebug)
        {
            obsVectors.Clear();
            obsVectors.Append("\nChild Dist Direct\n");
            obsVectors.Append("\nChild Pos\n");
            obsVectors.Append(layoutGrid.Centers[childA]);
            obsVectors.Append(layoutGrid.Centers[childB]);
            obsVectors.Append("\ndirection\n");
            obsVectors.Append(direction);
            obsVectors.Append("\nChild Index\n");
            obsVectors.Append(childA + " " + childB);

            obsVectors.Append("\n");
            obsVectors.Append("\nDirection Obs\n");
            obsVectors.Append(rightObs[0]);
            obsVectors.Append(" ");
            obsVectors.Append(rightObs[1]);
            obsVectors.Append("\n");
            obsVectors.Append(leftObs[0]);
            obsVectors.Append(" ");
            obsVectors.Append(leftObs[1]);
            obsVectors.Append("\n");
            obsVectors.Append(upObs[0]);
            obsVectors.Append(" ");
            obsVectors.Append(upObs[1]);
            obsVectors.Append("\n");
            obsVectors.Append(downObs[0]);
            obsVectors.Append(" ");
            obsVectors.Append(downObs[1]);

            obsVectors.Append("\n");
            obsVectors.Append("\nRatioDiff\n");
            obsVectors.Append(ratioDiff.ToString("F3"));
            obsVectors.Append("\nreward\n");
            float reward = CalucReward();
            obsVectors.Append(reward.ToString("F3"));

            obsVectors.Append("\n");
            obsVectors.Append("\nlayoutGrid.targetRatio\n");
            obsVectors.Append(layoutGrid.targetRatio[agentNum].ToString("F3"));
            obsVectors.Append("\n");
            obsVectors.Append(layoutGrid.currentRatio[agentNum].ToString("F3"));


        }

        SetMask();
    }

    private float CalucReward()
    {
        layoutGrid.MakeCurrentRatio();
        float ratioDiff = layoutGrid.targetRatio[agentNum] - layoutGrid.currentRatio[agentNum];//この差分が小さいほど評価が高い

        float targetRatio = layoutGrid.targetRatio[agentNum];//ratio is 0 to one so, ratioDiff range is 0 to targetRatio

        float xVal = (float)Remap((Mathf.Abs(targetRatio - Mathf.Abs(ratioDiff))), 0, targetRatio, -1, 1);
        float yVal = Mathf.Tan(xVal);
     
        return (float)Remap(yVal, -1.55, 1.55, 0, 1);
        //ロジット関数 NaNがでる
        /*xVal = Mathf.Log(xVal) - Mathf.Log(1f - xVal);

        if (xVal > 5f)
        {
            xVal = 5f;
        }
        if (xVal < -5f)
        {
            xVal = -5f;
        }

        //xVal = Mathf.Clamp(xVal, -5f, 5f) + 5f; // 0-10
        xVal += 0.5f;
        xVal *= 0.1f; // 0-1
       // xVal = (float)Remap(xVal, 0f, 10.0f, 0f, 1.0f);
         */

        //シグモイド関数
        // float xVal = (float)Remap((Mathf.Abs(targetRatio - Mathf.Abs(ratioDiff))), 0, targetRatio, -1, 1);
        // float reward = 1f / (1f + Mathf.Pow(2.718f, xVal * -4));


    }

    public static double Remap(double value, double from1, double to1, double from2, double to2)
    {
        if ((value - from1) == 0 && (to1 - from1) == 0)
        {
            return from2;
        }
        return (value - from1) / (to1 - from1) * (to2 - from2) + from2;
    }
    /// <summary>
    /// Status が Deadナラ　ソノ の方向にMaskをする
    /// </summary>
    private void SetMask()
    {
        List<Axis> A_OpenNeighbor = layoutGrid.GetOpenNeighbor(layoutGrid.AgentChild[agentNum][(int)Child.A]);

        if (!A_OpenNeighbor.Contains(Axis.Right))
        {
            SetActionMask((int)Child.A , (int)Axis.Right);
        }

        if (!A_OpenNeighbor.Contains(Axis.Left))
        {
            SetActionMask((int)Child.A, (int)Axis.Left);
        }

        if (!A_OpenNeighbor.Contains(Axis.Up))
        {
            SetActionMask((int)Child.A, (int)Axis.Up);
        }

        if (!A_OpenNeighbor.Contains(Axis.Down))
        {
            SetActionMask((int)Child.A, (int)Axis.Down);
        }


        List<Axis> B_OpenNeighbor = layoutGrid.GetOpenNeighbor(layoutGrid.AgentChild[agentNum][(int)Child.B]);

        if (!B_OpenNeighbor.Contains(Axis.Right))
        {
            SetActionMask((int)Child.B, (int)Axis.Right);
        }

        if (!B_OpenNeighbor.Contains(Axis.Left))
        {
            SetActionMask((int)Child.B, (int)Axis.Left);
        }

        if (!B_OpenNeighbor.Contains(Axis.Up))
        {
            SetActionMask((int)Child.B, (int)Axis.Up);
        }

        if (!B_OpenNeighbor.Contains(Axis.Down))
        {
            SetActionMask((int)Child.B, (int)Axis.Down);
        }
    }
   

    public override void AgentAction(float[] vectorAction)
    {
        Axis axisA = Mathf.FloorToInt(vectorAction[0]).ToAxis<int>();
        layoutGrid.AgentMove(agentNum, Child.A, axisA);

        Axis axisB = Mathf.FloorToInt(vectorAction[1]).ToAxis<int>();
        layoutGrid.AgentMove(agentNum, Child.B, axisB);

        SetReward(CalucReward());
    }

    public override void AgentReset()
    {
        if (transform.name == "0")
        {
            layoutEnviroment.ResetEnviroment(showDebug);
        }
    }

    public override float[] Heuristic()
    {
        var rtnList = new float[2];

        List<Axis> openNeighbor = layoutGrid.GetOpenNeighbor(layoutGrid.AgentChild[agentNum][(int)Child.A]);
        if (openNeighbor.Count == 0)
        {
            rtnList[(int)Child.A] = (int)Axis.None;
        }
        else
        {
            openNeighbor.Jitter();
            rtnList[(int)Child.A] = (int)openNeighbor[0];
        }

        openNeighbor = layoutGrid.GetOpenNeighbor(layoutGrid.AgentChild[agentNum][(int)Child.B]);
        if (openNeighbor.Count == 0)
        {
            rtnList[(int)Child.B] = (int)Axis.None;
        }
        else
        {
            openNeighbor.Jitter();
            rtnList[(int)Child.B] = (int)openNeighbor[0];
        }
        
        return rtnList;
    }
  
    private void OnGUI()
    {
        if (showDebug)
        {
           // CollectObservations();
            layoutGrid.PreviewOwner(layoutEnviroment.agentColors);
            myGUIUtilty.TextItemize(new Vector3(-5, layoutGrid.Y_Extend), obsVectors, Color.black);
        }

        if (!Academy.Instance.IsCommunicatorOn)
        {
            float targetRatio = layoutGrid.targetRatio[agentNum];
            float currentRatio = layoutGrid.currentRatio[agentNum];
            int subGridNum = Mathf.FloorToInt(layoutGrid.AliveIndex.Length * targetRatio) - Mathf.FloorToInt(layoutGrid.AliveIndex.Length * currentRatio);
         

            var centerA = layoutGrid.Centers[layoutGrid.AgentChild[agentNum][(int)Child.A]];
            var centerB = layoutGrid.Centers[layoutGrid.AgentChild[agentNum][(int)Child.B]];
            var avePos = (centerA + centerB) / 2f;
            
            myGUIUtilty.TextItemize(avePos, subGridNum.ToString(), Color.black);
        }
    }
}
