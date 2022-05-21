import streamlit as st
import pandas as pd
import numpy as np
from Bio import Align
from Bio import SeqIO
from pandas import DataFrame
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    st.set_page_config(layout="wide")
    st.session_state['expectedFolder'] = os.path.join('data')
    st.session_state['distanceData'] = os.path.join('data', 'fastq')
    if ('expectedFolder' not in st.session_state or 'distanceData' not in st.session_state
            or not os.path.exists(st.session_state['expectedFolder']) or not os.path.exists(
                st.session_state['distanceData'])):
        with st.expander("Data input paths"):
            st.session_state['expectedFolder'] = st.text_input("Expected folder data: ")
            st.session_state['distanceData'] = st.text_input("Test folder data: ")


    else:
        expectedFolder = {}
        dataFolder = {}
        algorithmType = ['local','global']

        aligner = Align.PairwiseAligner()

        for file in os.listdir(st.session_state['expectedFolder']):
            if (os.path.isfile(os.path.join(st.session_state['expectedFolder'], file))):
                expectedFolder[file.split('_')[0]] = os.path.join(st.session_state['expectedFolder'], file)
        with st.sidebar.expander('Expected file data: '):
            exFile = st.radio("", expectedFolder.keys(), index=0)

        for file in os.listdir(st.session_state['distanceData']):
            if (os.path.isfile(os.path.join(st.session_state['distanceData'], file))):
                dataFolder[file.split('_')[0]] = os.path.join(st.session_state['distanceData'], file)
        with st.sidebar.expander('File to calculate distance: '):
            dataFile = st.radio("", dataFolder.keys(), index=18)

        with st.sidebar.expander('Which algorithm to use: '):
            algorithmTy = st.radio("", algorithmType, index=1)

        with st.sidebar.expander('Aligner scores: '):
            aligner.match_score = st.number_input('Match score: ', value=1.0)
            aligner.mismatch_score = st.number_input('Mismatch score: ', value=0.0)
            aligner.open_gap_score = st.number_input('Open gap score: ', value=-1.0)
            aligner.extend_gap_score = st.number_input('Extend gap score: ', value=-1.0)
            aligner.target_end_gap_score = st.number_input('Target end gap score: ', value=0.0)
            aligner.query_end_gap_score = st.number_input('Query end gap score: ', value=-0.5)

        expected = list(SeqIO.parse(expectedFolder[exFile], "fasta"))
        data = list(SeqIO.parse(dataFolder[dataFile], "fastq"))

        cols = st.columns(len(expected))
        distances = {}
        for i in range(len(expected)):
            distances[i] = []
            with cols[i]:
                with st.expander('Sequence ' + expected[i].id):
                    exp = expected[i].seq
                    for sequence in data:
                        alignments = aligner.align(expected[i].seq, sequence.seq)
                        distance = int(min(len(exp), len(sequence.seq)) - alignments[0].score)
                        distances[i].append(distance)
                    fig, ax = plt.subplots()
                    ax.hist(distances[i], bins=30)
                    st.pyplot(fig)

                    fig, ax = plt.subplots()
                    ax.boxplot(distances[i])
                    st.pyplot(fig)