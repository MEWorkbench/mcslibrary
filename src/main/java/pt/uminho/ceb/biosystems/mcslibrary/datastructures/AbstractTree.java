package pt.uminho.ceb.biosystems.mcslibrary.datastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class AbstractTree<T extends Object> {
	private List<AbstractTree<T>> children;
	private T nodeValue;
	private AbstractTree<T> parent;
	
	public AbstractTree() {
		this.parent = null;
		this.children = new ArrayList<>();
		this.nodeValue = null;
	}
	
	public AbstractTree(T value) {
		this.children = new ArrayList<>();
		this.nodeValue = value;
		this.parent = null;
	}
	
	public AbstractTree(T value, AbstractTree<T> parent) {
		this.children = new ArrayList<>();
		this.nodeValue = value;
		this.parent = parent;
	}
	
	public boolean isLeaf(){
		return getNumberOfChildren() == 0;
	}
	
	public T getValue(){
		return this.nodeValue;
	}
	
	public AbstractTree<T> getParent(){
		return this.parent;
	}
	
	public void setValue(T value){
		this.nodeValue = value;
	}
	
	public void setParent(AbstractTree<T> parent){
		this.parent = parent;
	}
	
	public AbstractTree<T> getChild(int index){
		return this.children.get(index);
	}
	
	public int getNumberOfChildren(){
		return this.children.size();
	}
	
	public void addChild(T value){
		this.children.add(new AbstractTree<T>(value, this));
	}
	
	public void addChild(AbstractTree<T> node){
		node.parent = this;
		this.children.add(node);
	}
	
	public void removeChild(int index){
		this.children.remove(index);
	}
	
	public void removeChild(AbstractTree<T> node){
		this.children.remove(node);
	}
	
	public List<AbstractTree<T>> getAllNodes(){
		ArrayList<AbstractTree<T>> l = new ArrayList<AbstractTree<T>>();
		getAllNodes(l);
		return l;
	}
	
	public void getAllNodes(List<AbstractTree<T>> list){
		list.addAll(children);
		for (int i = 0; i < getNumberOfChildren(); i++) {
			children.get(i).getAllNodes(list);
		}
	}
	
	public void removeAllChildren(Collection<AbstractTree<T>> children){
		this.children.removeAll(children);
	}
	
	public String toString(){
		return printRecursive();
	}
	
	private String printRecursive(){
		return getValue()+"\n"+printChildren(1);
	}
	
	public List<AbstractTree<T>> getPath(){
		List<AbstractTree<T>> list = new ArrayList<AbstractTree<T>>();
		getPath(list);
		return list;
	}
	
	
	private void getPath(List<AbstractTree<T>> list){
//		List<AbstractTree<T>> path = new ArrayList<AbstractTree<T>>();
		if (this.parent != null) {
			this.parent.getPath(list);
			list.add(this);
		}

		
	}
	private String printChildren(int curLevel){
		String str = "";
		for (int i = 0; i < getNumberOfChildren(); i++) {
			str += repeat(curLevel, "\t")+getChild(i).getValue()+"\n";
			str += getChild(i).printChildren(curLevel+1);
		}
		return str;
	}
	
	private static String repeat(int count, String with) {
	    return new String(new char[count]).replace("\0", with);
	}

	private static String repeat(int count) {
	    return repeat(count, " ");
	}

	public List<AbstractTree<T>> getAllLeaves(){
		List<AbstractTree<T>> list = new ArrayList<AbstractTree<T>>();
		getAllLeaves(list);
		return list;
	}
	
	private void getAllLeaves(List<AbstractTree<T>> list){
		if (!isLeaf()) {
			for (int i = 0; i < getNumberOfChildren(); i++) {
				getChild(i).getAllLeaves(list);
			}
		} else {
			list.add(this);
		}
	}
}
