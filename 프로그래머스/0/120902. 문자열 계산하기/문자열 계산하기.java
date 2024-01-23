class Solution {
    public int solution(String my_string) {
        String[] sp = my_string.split("\\s");
        int a = Integer.parseInt(sp[0]);
        for(int i = 1; i < sp.length; i += 2){
            int b = Integer.parseInt(sp[i + 1]);
            a += sp[i].equals("+") ? b : -b;
        }
        return a;
    }
}