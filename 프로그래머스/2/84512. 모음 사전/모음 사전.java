class Solution {
    char[] aeiou = new char[] { 'A', 'E', 'I', 'O', 'U' };
    int count = 0;
    String word;
    boolean search = false;
    public int solution(String word) {
        this.word = word;
        dfs("");
        return count;
    }
    public void dfs(String s){
        if(search || s.length() > 5) return;
        else if(word.equals(s)){
            search = true;
            return;
        }else{
            count++;
            for(int i = 0; i < aeiou.length; i++) dfs(s + aeiou[i]);
        }
    }
}