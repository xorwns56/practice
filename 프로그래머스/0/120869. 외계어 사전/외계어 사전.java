class Solution {
    public int solution(String[] spell, String[] dic) {
        int answer = 0;
        for(int i = 0; i < dic.length; i++){
            String s = dic[i];
            boolean contains = true;
            for(int j = 0; j < spell.length; j++){
                if(!s.contains(spell[j])){
                    contains = false;
                    break;
                }
                s = s.replaceFirst(spell[j], "");
            }
            if(contains && s.isEmpty()) return 1;
        }
        return 2;
    }
}