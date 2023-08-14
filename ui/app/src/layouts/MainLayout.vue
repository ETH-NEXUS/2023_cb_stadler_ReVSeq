<script setup lang="ts">
import {ref} from 'vue'
import {useUserStore} from 'stores/user'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'

const router = useRouter()

const sideBarOpen = ref<boolean>(false)
const userStore = useUserStore()
const {t} = useI18n()

const toggleSideBar = () => {
  sideBarOpen.value = !sideBarOpen.value
}

const logoutButtonClicked = async () => {
  try {
    await userStore.sessionLogout()
    await router.push({path: '/login'})
  } catch (err) {
    console.error(err)
  }
}

const navigateToSearch = async () => {
  try {
    await router.push({path: '/search'})
  } catch (err) {
    console.error(err)
  }
}

const navigateToAboutPage = async () => {
  try {
    await router.push({path: '/about'})
  } catch (err) {
    console.error(err)
  }
}
const navigateToHome = async () => {
  try {
    await router.push({path: '/'})
  } catch (err) {
    console.error(err)
  }
}

const links = [
  {name: 'Home', function: () => navigateToHome()},
  {name: 'Database', function: () => navigateToSearch()},
  {name: 'About', function: () => navigateToAboutPage()},
]
</script>

<template>
  <q-layout view="lHh Lpr lFf" class="tw-font-sans">
    <section class="tw-relative">
      <div class="tw-absolute tw-top-0 tw-right-0 tw-flex tw-w-full tw-h-3/4 md:tw-h-2/3">
        <img
          class="tw-w-64 md:tw-w-80 2xl:tw-w-auto tw-self-start"
          src="../assets/pattern-dots-indigo-left.svg"
          alt="orange background 1" />
        <img
          class="tw-w-64 md:tw-w-80 2xl:tw-w-auto tw-self-end tw-ml-auto"
          src="../assets/pattern-dots-indigo-right.svg"
          alt="orange background 2" />
      </div>
      <nav class="tw-flex tw-mb-20 tw-justify-between tw-items-center tw-py-6 tw-px-10 tw-relative">
        <a class="tw-text-lg tw-font-bold">
          <img class="tw-h-8" src="../assets/quasar-logo-vertical.svg" alt="logo image" width="auto" />
        </a>
        <div class="xl:tw-hidden">
          <button
            @click="toggleSideBar"
            class="tw-navbar-burger focus:tw-outline-none tw-text-indigo-900 hover:tw-text-indigo-800">
            <svg
              class="tw-block tw-h-6 tw-w-6"
              fill="currentColor"
              viewbox="0 0 20 20"
              xmlns="http://www.w3.org/2000/svg">
              <title>Mobile menu</title>
              <path d="M0 3h20v2H0V3zm0 6h20v2H0V9zm0 6h20v2H0v-2z"></path>
            </svg>
          </button>
        </div>
        <ul
          class="tw-hidden xl:tw-flex tw-absolute tw-top-1/2 tw-left-1/2 tw-transform tw--translate-x-1/2 tw--translate-y-1/2">
          <li v-for="(link, index) in links" :key="link.name + index">
            <a
              class="cursor-pointer tw-text-lg tw-mr-10 2xl:tw-mr-16 tw-font-extrabold hover:tw-text-indigo-800"
              @click="link.function">
              {{ link.name }}
            </a>
          </li>
        </ul>
        <div class="tw-hidden xl:tw-flex tw-items-center">
          <a
            @click="logoutButtonClicked"
            class="cursor-pointer tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200">
            {{ t('label.logout') }}
          </a>
        </div>
      </nav>
      <div class="tw-container tw-px-4 tw-mx-auto tw-relative">
        <div class="tw-max-w-5xl tw-mx-auto tw-text-center">
          <span class="tw-text-xl md:tw-text-2xl tw-font-extrabold tw-text-orange-500">
            {{ t('main_layout.eth') }}
          </span>

          <h2
            class="tw-max-w-4xl tw-mx-auto tw-text-3xl sm:tw-text-4xl lg:tw-text-5xl tw-font-extrabold tw-font-heading tw-mt-1 tw-mb-6">
            {{ t('main_layout.title') }}
          </h2>
          <p class="tw-text-xl md:tw-text-2xl tw-font-extrabold tw-leading-8 tw-mb-12">
            {{ t('main_layout.short_description') }}
          </p>
          <!--          <p class="tw-text-xl md:tw-text-2xl tw-font-extrabold tw-leading-8 tw-mb-12">-->
          <!--            {{ t('main_layout.description') }}-->
          <!--          </p>-->
          <div class="tw-flex tw-flex-wrap tw-mb-20 tw-justify-center">
            <a
              @click="navigateToSearch"
              class="cursor-pointer tw-inline-block tw-w-full md:tw-w-auto tw-mb-2 md:tw-mb-0 md:tw-mr-4 tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200">
              {{ t('label.search_db') }}
            </a>
            <a
              @click="navigateToAboutPage"
              class="cursor-pointer tw-inline-block tw-w-full md:tw-w-auto tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg hover:tw-text-white tw-font-extrabold tw-bg-white tw-border-3 tw-border-indigo-900 hover:tw-bg-indigo-800 tw-shadow tw-rounded tw-transition tw-duration-200">
              {{ t('label.about') }}
            </a>
          </div>
        </div>
      </div>

      <div class="tw-navbar-menu tw-relative tw-z-50" :class="sideBarOpen ? '' : 'tw-hidden'">
        <div class="tw-navbar-backdrop tw-fixed tw-inset-0 tw-bg-gray-800 tw-opacity-25"></div>
        <nav
          class="tw-fixed tw-top-0 tw-left-0 tw-bottom-0 tw-flex tw-flex-col tw-w-full md:tw-w-5/6 tw-max-w-sm tw-py-8 tw-px-8 tw-bg-white tw-border-r tw-overflow-y-auto">
          <div class="tw-flex tw-items-center tw-mb-8">
            <a class="tw-mr-auto tw-text-2xl tw-font-bold tw-leading-none">
              <img class="tw-h-6" src="../assets/quasar-logo-vertical.svg" alt="" width="auto" />
            </a>
            <button class="tw-navbar-close" @click="toggleSideBar">
              <svg
                class="tw-h-6 tw-w-6 tw-text-gray-500 tw-cursor-pointer hover:tw-text-gray-500"
                xmlns="http://www.w3.org/2000/svg"
                fill="none"
                viewbox="0 0 24 24"
                stroke="currentColor">
                <path
                  stroke-linecap="round"
                  stroke-linejoin="round"
                  stroke-width="2"
                  d="M6 18L18 6M6 6l12 12"></path>
              </svg>
            </button>
          </div>
          <div class="tw-my-auto">
            <ul class="tw-py-10">
              <li class="tw-mb-1" v-for="link in links" :key="link.name">
                <a
                  class="cursor-pointer tw-block tw-p-4 tw-text-lg tw-font-extrabold hover:tw-bg-gray-50 rounded"
                  @click="link.function">
                  {{ link.name }}
                </a>
              </li>
            </ul>
          </div>
          <div>
            <a
              @click="logoutButtonClicked"
              class="tw-block tw-py-4 tw-px-6 ttw-ext-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200">
              {{ t('label.logout') }}
            </a>
          </div>
        </nav>
      </div>
      <q-page-container>
        <router-view />
      </q-page-container>
    </section>
    <section>
      <div class="tw-mt-36 tw-mb-8">
        <div class="tw-container tw-px-4 tw-mx-auto tw-text-center">
          <a class="tw-inline-block tw-mx-auto tw-mb-2" href="https://www.nexus.ethz.ch/events.html">
            <img class="tw-h-12" src="../assets/nexus_logo.png" alt="" />
          </a>
        </div>
        <div class="tw-container tw-px-4 tw-mx-auto">
          <p class="tw-text-center tw-text-md md:tw-text-xl tw-font-extrabold">
            {{ `${new Date().getFullYear()} ETH Zurich` }}
          </p>
        </div>
      </div>
    </section>
  </q-layout>
</template>

<style>
@import 'tailwindcss/base.css';
@import 'tailwindcss/components.css';
@import 'tailwindcss/utilities.css';
</style>
