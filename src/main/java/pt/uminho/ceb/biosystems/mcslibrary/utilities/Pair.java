/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
package pt.uminho.ceb.biosystems.mcslibrary.utilities;

public class Pair<T1, T2> {

	private T1 a;
	private T2 b;

	public Pair(T1 a, T2 b) {
		this.a = a;
		this.b = b;
	}

	public T1 getA() {
		return a;
	}

	public T2 getB() {
		return b;
	}

	public void setA(T1 a) {
		this.a = a;
	}

	public void setB(T2 b) {
		this.b = b;
	}
	
	public String toString(){
		String str = "("+a+"; "+b+")";
		return str;
	}
}
