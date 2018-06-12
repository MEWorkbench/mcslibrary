package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RegexUtils {
	
	public static List<List<String>> findAll(String pattern, String toMatch){
		List<List<String>> res = new ArrayList<List<String>>();
		Pattern p = Pattern.compile(pattern);
		Matcher m = p.matcher(toMatch);
		while (m.find()){
			List<String> cur = new ArrayList<String>();
			for (int i = 0; i < m.groupCount(); i++) {
				cur.add(m.group(i));
			}
			res.add(cur);
		}
		return res;
	}
}
