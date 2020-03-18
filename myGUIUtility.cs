using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Text;
public  class myGUIUtilty : MonoBehaviour
{
    public static void TextItemize(Rect rect, List<string> texts, Color color, string title = "", int fontSize = 10)
    {

        StringBuilder mergeText = new StringBuilder();
        if (title != "")
        {
            mergeText.Append(title);
            mergeText.Append("\n");
        }
        foreach (var text in texts)
        {
            mergeText.Append(text);
            mergeText.Append("\n");
        }

        var gUIContent = new GUIContent(mergeText.ToString());

        var guiStyle = new GUIStyle();
        guiStyle.fontSize = fontSize;
        guiStyle.CalcSize(gUIContent);
        guiStyle.normal.textColor = color;


        GUI.TextArea(rect, mergeText.ToString(), guiStyle);
    }

 
    public static void TextItemize(Vector3 pos, List<string> texts, Color color, string title = "", int fontSize = 10)
    {
        pos = Camera.main.WorldToScreenPoint(pos);
        pos.y = Screen.height - pos.y;
        var rect = new Rect(pos, Vector2.one);
     
        StringBuilder mergeText = new StringBuilder();
        if (title != "")
        {
            mergeText.Append(title);
            mergeText.Append("\n");
        }
        foreach (var text in texts)
        {
            mergeText.Append(text);
            mergeText.Append("\n");
        }

        var gUIContent = new GUIContent(mergeText.ToString());

        var guiStyle = new GUIStyle();
        guiStyle.fontSize = fontSize;
        guiStyle.CalcSize(gUIContent);
        guiStyle.normal.textColor = color;

        GUI.TextArea(rect, mergeText.ToString(), guiStyle);
    }

    public static void TextItemize(Vector3 pos, StringBuilder text, Color color, int fontSize = 10)
    {
        pos = Camera.main.WorldToScreenPoint(pos);
        pos.y = Screen.height - pos.y;
        var rect = new Rect(pos, Vector2.one);

        var gUIContent = new GUIContent(text.ToString());

        var guiStyle = new GUIStyle();
        guiStyle.fontSize = fontSize;
        guiStyle.CalcSize(gUIContent);
        guiStyle.normal.textColor = color;

        GUI.TextArea(rect, text.ToString(), guiStyle);
    }
    public static void TextItemize(Vector2 pos, StringBuilder text, Color color, int fontSize = 10)
    {
        pos = Camera.main.WorldToScreenPoint(pos);
        pos.y = Screen.height - pos.y;
        var rect = new Rect(pos, Vector2.one);

        var gUIContent = new GUIContent(text.ToString());

        var guiStyle = new GUIStyle();
        guiStyle.fontSize = fontSize;
        guiStyle.CalcSize(gUIContent);
        guiStyle.normal.textColor = color;

        GUI.TextArea(rect, text.ToString(), guiStyle);
    }

    public static void TextItemize(Vector2 pos, string text, Color color, int fontSize = 10)
    {
        pos = Camera.main.WorldToScreenPoint(pos);
        pos.y = Screen.height - pos.y;
        var rect = new Rect(pos, Vector2.one);

        var gUIContent = new GUIContent(text);

        var guiStyle = new GUIStyle();
        guiStyle.fontSize = fontSize;
        guiStyle.CalcSize(gUIContent);
        guiStyle.normal.textColor = color;

        GUI.TextArea(rect, text, guiStyle);
    }
}