class Solution {
    public String solution(String letter) {
        String[] morse = {".-","-...","-.-.","-..",".","..-.","--.","....","..",".---","-.-",".-..","--","-.","---",".--.","--.-",".-.","...","-","..-","...-",".--","-..-","-.--","--.."}; 
        String[] sp = letter.split("\\s");
        String answer = "";
        for(int i = 0; i < sp.length; i++){
            for(int j = 0; j < morse.length; j++){
                if(sp[i].equals(morse[j])){
                    answer += (char)(j + 'a');
                    break;
                }
            }
        }
        return answer;
    }
}