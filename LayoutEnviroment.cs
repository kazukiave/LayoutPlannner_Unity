using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class LayoutEnviroment : MonoBehaviour
{

    GameObject[] agents;
    [System.NonSerialized]
    public Color[] agentColors;
    [System.NonSerialized]
    public LayoutGrid layoutGrid;

    public int reduceNum = 0;
    public int offsetNum = 1;
    public float scatterRate = 1.1f;
    public float ratioSum = 1.0f;

    public int gridSize;
    public int xExtend;
    public int yExtend;

    [System.NonSerialized]
    public bool Init = false;

    public void InitEnviroment()
    {
        if (Init) return;

        //子のエージェントを参照
        var transforms = GetComponentsInChildren<Transform>(false);
        transforms = transforms.Where(trans => trans.CompareTag("agent")).ToArray();
        agents = new GameObject[transforms.Length];
        for (int i = 0; i < transforms.Length; i++)
        {
            transforms[i].gameObject.name = i.ToString();
            agents[i] = transforms[i].gameObject;
            agents[i].transform.name = i.ToString();
        }

        //make colors for Preview
        agentColors = new Color[agents.Length];
        for (int i = 0; i < agents.Length; i++)
        {
            agentColors[i] = Color.HSVToRGB(((float)i / (float)(agents.Length)), 1f, 1f, true);
        }


        //layoutGridを初期化
        layoutGrid = gameObject.AddComponent<LayoutGrid>();
        layoutGrid = new LayoutGrid(gridSize, xExtend, yExtend, agents.Length);

        ResetEnviroment();

        Init = true;
    }

    
    public void ResetEnviroment(bool ShowDebug = false)
    {
        layoutGrid.ResetGrid();
      
        //Gridの形状を決める
        layoutGrid.RandomReduce(reduceNum, offsetNum);
        if(ShowDebug) layoutGrid.PreviewStatues();
        layoutGrid.MakeTargetRatio(scatterRate, ratioSum);
       
        //Agentの初期位置を決める。
        var firstPos = layoutGrid.FirstPosition();
        layoutGrid.SplitPosition(firstPos);
    }


}
