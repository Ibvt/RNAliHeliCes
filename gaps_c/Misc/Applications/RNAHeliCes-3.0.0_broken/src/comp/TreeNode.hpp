//*****************************************************************************
#ifndef __TREENODE_H
#define __TREENODE_H 1

/**
@file TreeNode.h
****************
@brief TreeNode.h is the header file for @c val::biocpp::tools::CTreeNode
class.

TreeNode.h is the header file for a basic tree node template
@c val::biocpp::tools::CTreeNode. Tree nodes can be used in the case of tree
comparison. A tree node can have at most one parent and as many children as
needed.

@see val::biocpp::tools::CTreeNode

@author Valentin GUIGNON
@version 1.0
@date 06/09/2004
*/
#include <list>
#include "DebugTools.hpp"

namespace val
{
namespace biocpp
{
namespace tools
{
/**
@class CTreeNode
****************
@brief A basic tree node template.

This templates describes the behaviour of a tree node containing an object.
Each node creates its own copy of the object.

@note the object can be a pointer to something managed by the programmer.
Parent nodes manage memory allocated for child nodes (ie.: the programmer can
create nodes but only need to delete the parent node).

@note this class prevents cyclic-trees: a parent node can't be added to
any of its children or sub*children.

@b Templates: \n
@param Object: \n
 tree object type.\n
 Constraints: the type must support the equal operator. It can be a pointer but
  nor an array.

@author Valentin GUIGNON
@version 1.3
@date 21/08/2004
*/
template <class Object>
class CTreeNode
{
// public typdefs
public:
	//! Tree node pointer type.
	typedef CTreeNode<Object>*				CTreeNodePointer;
	//! Constant tree node pointer type.
	typedef const CTreeNode<Object>*		ConstCTreeNodePointer;
	//! List of tree node pointers type.
	typedef std::list<CTreeNodePointer>		CTreeNodePointerList;


// protected members
protected:
	//! object owned by the node
	Object								m_tObject;
	//! pointer to a parent node if one (NULL otherwise)
	CTreeNodePointer					m_ptParent;
	//! children nodes lists
	CTreeNodePointerList				m_tChildrenList;


// protected methods
protected:
	/**
	@fn void FillSuffixList(CTreeNodePointerList& tNodesList)
	*********************************************************
	@brief Fills a suffix list recursively.

	It pushes in the end of the given list, the children encountered using the
	suffix order in the sub-tree starting from this node.

	@param CTreeNodePointerList& tNodesList:\n
	 the list to fill.
	*/
	void FillSuffixList(CTreeNodePointerList& tNodesList)
	{
		TRACEI("void FillSuffixList(CTreeNodePointerList& tNodesList)");
		// follows the children nodes
		for (typename CTreeNodePointerList::const_iterator tIter = m_tChildrenList.begin(); tIter != m_tChildrenList.end(); tIter++)
		{
			(*tIter)->FillSuffixList(tNodesList);
		}
		// adds this node to the list
		tNodesList.push_back(this);
	}


public:
	/**
	@fn CTreeNode(const Object &tObject)
	************************************
	@brief Standard constructor.

	Standard constructor.

	@param const Object &tObject:\n
	 object for this node.
	*/
	CTreeNode(const Object &tObject)
	: m_tObject(tObject)
	, m_ptParent(NULL)
	{
		TRACEI("CTreeNode(const Object &tObject)");
	}


	/**
	@fn ~CTreeNode()
	****************
	@brief Standard virtual destructor.

	Standard virtual destructor.
	*/
	virtual ~CTreeNode()
	{
		TRACEI("~CTreeNode()");
		// remove from parent
		if (NULL != m_ptParent)
		{
			m_ptParent->m_tChildrenList.remove(this);
		}
		// delete children list
		for (typename CTreeNodePointerList::iterator tIter = m_tChildrenList.begin(); tIter != m_tChildrenList.end(); tIter++)
		{
			(*tIter)->m_ptParent = NULL;
			delete (*tIter);
		}
	}


	/**
	@fn const Object& GetObject() const
	***********************************
	@brief Returns node object.

	Returns node object.

	@return const Object&:\n
	 Object value for this node.
	*/
	const Object& GetObject() const
	{
		TRACEI("const Object& GetObject() const");
		return m_tObject;
	}


	/**
	@fn bool HasParent() const
	**************************
	@brief Tells if the node has a parent.

	Returns true if this node has a parent.

	@return bool:\n
	 true is the node has a parent, false otherwise.
	*/
	bool HasParent() const
	{
		TRACEI("bool HasParent() const");
		return (NULL == m_ptParent)? false:true;
	}


	/**
	@fn bool HasChildren() const
	****************************
	@brief Tells if the node has children (or at least one child).

	Returns true if this node has at least a child.

	@return bool:\n
	 true is the node has at least a child, false otherwise.
	*/
	bool HasChildren() const
	{
		TRACEI("bool HasChildren() const");
		return m_tChildrenList.empty()? false:true;
	}


	/**
	@fn bool IsDescendantOf(const CTreeNodePointer &ptTreeNode) const
	*****************************************************************
	@brief Tells if current node is the decendant of the given node.

	Returns true if current node is a child (or a subchild) of the given node.

	@param const CTreeNodePointer &ptTreeNode:\n
	 parent node to check.

	@return bool:\n
	 returns true if current node is a descendant of the given node.
	*/
	bool IsDescendantOf(const CTreeNodePointer &ptTreeNode) const
	{
		TRACEI("bool IsDescendantOf(const CTreeNodePointer &ptTreeNode) const");
		bool bRet = false;
		if ((NULL != m_ptParent) && (NULL != ptTreeNode))
		{
			if (ptTreeNode == m_ptParent)
			{
				bRet = true;
			}
			else
			{
				bRet = m_ptParent->IsDescendantOf(ptTreeNode);
			}
		}
		return bRet;
	}


	/**
	@fn bool IsAncestorOf(const CTreeNodePointer &ptTreeNode) const
	***************************************************************
	@brief Tells if current node is the ancestor of the given node.

	Returns true if current node is an ancestor of the given node.

	@param const CTreeNodePointer &ptTreeNode:\n
	 child node to check.

	@return bool:\n
	 returns true if current node is an ancestor of the given node.
	*/
	bool IsAncestorOf(const CTreeNodePointer &ptTreeNode) const
	{
		TRACED("bool IsAncestorOf(const CTreeNodePointer &ptTreeNode) const");
		bool bRet = false;
		if (NULL != ptTreeNode)
		{
			bRet = ptTreeNode->IsDescendantOf(this);
		}
		return bRet;
	}


	/**
	@fn bool IsDirectLeftChildOf(const CTreeNodePointer &ptTreeNode) const
	**********************************************************************
	@brief Tells if current node is the direct left child of the given node.

	Returns true if current node is the direct left child of the given node.

	@param const CTreeNodePointer &ptTreeNode:\n
	 parent node to check.

	@return bool:\n
	 returns true if current node is the direct left descendant of the given
	 node.
	*/
	bool IsDirectLeftChildOf(const CTreeNodePointer &ptTreeNode) const
	{
		TRACEI("bool IsDirectLeftChildOf(const CTreeNodePointer &ptTreeNode) const");
		bool bRet = false;
		if ((NULL != m_ptParent)
			&& (ptTreeNode == m_ptParent)
			&& (this == ptTreeNode->GetDirectChildrenList().front()))
		{
			bRet = true;
		}
		return bRet;
	}


	/**
	@fn CTreeNodePointer AddChild(const Object &tObject)
	****************************************************
	@brief Adds a child node to this node.

	Creates a child node using object tObject and adds it on specified side.

	@param const Object &tObject:\n
	 an object (will be copied) for the child node.

	@return CTreeNodePointer:\n
	 the pointer to the new node.
	*/
	CTreeNodePointer AddChild(const Object &tObject)
	{
		TRACEI("void AddChild(const Object &tObject)");
		CTreeNodePointer ptNewNode = new CTreeNode<Object>(tObject);
		ptNewNode->m_ptParent = this;
		m_tChildrenList.push_back(ptNewNode);
		return ptNewNode;
	}


	/**
	@fn int GetDirectChildrenCount() const
	**************************************
	@brief Gets the number of direct children of this node.

	Returns the count of direct children of this parent. This count does not
	include sub-children count.

	@return int:\n
	 number of children.
	*/
	int GetDirectChildrenCount() const
	{
		TRACEI("int GetDirectChildrenCount() const");
		return (signed)m_tChildrenList.size();
	}


	/**
	@fn int GetChildrenCount() const
	********************************
	@brief Gets the number of children and sub*children of this node.

	Returns the count of children and sub*children of this node.

	@return int:\n
	 number of children plus the children count of each children.
	*/
	int GetChildrenCount() const
	{
		TRACEI("int GetChildrenCount() const"); //+debug
		// adds direct child nodes count
		int iRet = (signed)m_tChildrenList.size();
		// get children count
		for (typename CTreeNodePointerList::const_iterator tIter = m_tChildrenList.begin(); tIter != m_tChildrenList.end(); tIter++)
		{
			iRet += (*tIter)->GetChildrenCount();
		}
		return iRet;
	}


	/**
	@fn CTreeNodePointer GetParentPointer()
	***************************************
	@brief Returns a pointer to parent node.

	Returns the parent node pointer if one, NULL otherwise.

	@return CTreeNodePointer:\n
	 parent node pointer or NULL.
	*/
	CTreeNodePointer GetParentPointer()
	{
		TRACEI("CTreeNodePointer GetParentPointer() const");
		return m_ptParent;
	}


	/**
	@fn CTreeNodePointer GetLeftLeafPointer()
	*****************************************
	@brief Returns a pointer to the lowest left child.

	Returns the lowest left child or sub-child if there are children, self
	otherwise.

	@return CTreeNodePointer:\n
	 the leftest leaf for this node or this node if there are no children.
	*/
	CTreeNodePointer GetLeftLeafPointer()
	{
		TRACEI("CTreeNodePointer GetLeftLeafPointer()");
		CTreeNodePointer ptLeftLeaf = this;
		if (0 < m_tChildrenList.size())
		{
			ptLeftLeaf = m_tChildrenList.front()->GetLeftLeafPointer();
		}
		return ptLeftLeaf;
	}


	/**
	@fn CTreeNodePointerList GetDirectChildrenList()
	************************************************
	@brief Returns the list of direct children pointers.

	Returns the pointer list of direct children. Direct children are children
	directly connected to this node.

	@return CTreeNodePointerList:\n
	 the list of direct children on the specified side.
	*/
	CTreeNodePointerList GetDirectChildrenList()
	{
		TRACEI("const CTreeNodePointerList& GetDirectChildrenList() const");
		return m_tChildrenList;
	}


	/**
	@fn CTreeNodePointerList GetSuffixList()
	****************************************
	@brief Returns the suffix list from this node.

	Returns a list containing all nodes below and this one ordered in a
	"suffix" order.

	@return CTreeNodePointerList:\n
	 an ordered list that contains all children and sub-children nodes
	 starting from this node. The front of the list contains the first node
	 encountered using the "suffix" way. Then each following node is the
	 following one using the suffix order.
	*/
	CTreeNodePointerList GetSuffixList()
	{
		TRACEI("const CTreeNodePointerList GetSuffixList() const");
		CTreeNodePointerList tNodesList;
		FillSuffixList(tNodesList);
		return tNodesList;
	}


/*
	Iterator getPrefixIterator()
	{
	}
	Iterator getInfixIterator()
	{
	}
*/
}; // class CTreeNode
}; // namespace tools
}; // namespace biocpp
}; // namespace val
#endif // #ifndef __TREENODE_H
